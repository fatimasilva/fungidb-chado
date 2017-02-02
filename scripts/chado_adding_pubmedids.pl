#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use DBI;
use IO::Prompt;

# Adds PubMed IDs to polypeptide features based on a list of GeneIDs. 
# Input file in tab format, with headers, first column gene IDs, second column PubMed IDs
# only one pair geneid:pubmedid per row
# Checks if the PubMed ID is already present
# For non codings genes (rRNAs and tRNAs) it adds the publication to the rRNA/tRNA feature 

my ($help, $f_data, $name_out, $dbname, $dbhost, $dbport, $dbuser);

GetOptions(
    'a|data_file=s'       => \$f_data,
    'o|output_prefix=s'   => \$name_out,
    'h|help=s'            => \$help,
    'd|database_name=s'   => \$dbname,
    's|database_host=s'   => \$dbhost,
    'u|database_user=s'   => \$dbuser,
    'p|database_port=s'   => \$dbport,

);

(($f_data && $name_out && $dbname && $dbuser && $dbhost && $dbport ) && !$help) || die <<USAGE;

Usage: $0
    -a|data_file <tab file with gene ID (first column) and pubmedid to load (second column)>
    -o|name_out <prefix for output files>
    -h|help    <This help message>
    -d|database_name
    -u|database_user
    -s|database_host
    -p|database_port
USAGE

#### log and output files
my $f_drop = $name_out . ".log";
my $f_out = $name_out . "_changes_chado.tsv";
my $f_temp = $name_out . ".temp";
my $f_tag = $name_out . "_pep_feature_ids_history.tsv";


#### Database connection
## Prompt user for db password
my $dbpass = prompt('Password:', -echo => '*');

$dbpass = "$dbpass";

## Connect to db
my $dbi_connect = "DBI:Pg:dbname=$dbname;host=$dbhost;port=$dbport";

my $dbh = DBI->connect($dbi_connect, $dbuser, $dbpass,{RaiseError => 0, AutoCommit => 0}) 
                        or die "Can't connect to db $dbname!\n"
                        ;

print "Success connecting to db $dbname!\n";

my $errflag;

#### Fungal instance
my $type_id_gene = 3596;
my $type_id_pseudo = 3229;
my $type_id_part = 8;
my $type_id_derives = 2852; 
my $type_id_unfetched = 2836;
my $type_id_rrna = 3145;
my $type_id_trna = 3146;
#my $type_id_sno = ; # not loaded yet, but they will

#### Pathogen instance
#my $type_id_part = 42;
#my $type_id_derives = 69; 
#my $type_id_unfetched = 26797;

#### Prepare SQL statements

## Get feature id of gene or pseudogene by gene ID (uses type_id_gene and type_id_pseudo
my $s_sql_get_gene_fid = $dbh->prepare('SELECT feature_id FROM feature WHERE uniquename=? AND (type_id=? OR type_id=?)');

## Get feature_id of mRNA and polypeptide (by feature relationships from gene ID)
my $s_sql_subject_id=$dbh->prepare('SELECT subject_id from feature_relationship where object_id=? and type_id=?');

## Get pub_id, the identifier of the PubMed ID in the pub table (uses type_id_unfetched)
my $s_sql_get_pub_id = $dbh->prepare('SELECT pub_id FROM pub WHERE lower(uniquename) = lower(?) AND type_id=?');

## Insert new pub id entry in the pub table (uses type_id_unfetched)
my $s_sql_add_pub_id = $dbh->prepare('INSERT INTO pub (uniquename, type_id, is_obsolete) VALUES (?,?,?)');

## Insert pub_id in feature_pub table
my $s_sql_add_feature_pub_id = $dbh->prepare('INSERT INTO feature_pub (feature_id, pub_id) VALUES (?,?)');

## Get feature_pub_id  
my $s_sql_get_fpub_id = $dbh->prepare('SELECT feature_pub_id FROM feature_pub WHERE feature_id=? AND pub_id=?');

## Check for non coding genes !00 only checks for rRNAs and tRNAs
my $s_sql_check_non_coding = $dbh->prepare('SELECT feature_id FROM feature WHERE feature_id=? AND type_id IN (?,?)');

#### Get PubMed IDs into hashes key:gene_id, value:pubmed_id 
open (my $h_drop, ">>", $f_drop);
open (my $h_temp, ">>", $f_temp);

my $pubmed_ids_ref = get_pubmed_ids_from_tsv($f_data);
my %pubmed_ids = %{$pubmed_ids_ref};

## Check extraction

foreach my $id ( sort ( keys %pubmed_ids) ){

    print $h_temp "$id\t$pubmed_ids{$id}\n";

}

open (my $h_out, ">>", $f_out);
open (my $h_tag, ">>" ,$f_tag);

#### Load pubmed_ids into chado 
Geneid: foreach my $gene_id ( sort (keys %pubmed_ids)){
    
    if (defined $errflag){
        warn "Error adding pubmed ID into feature: " . $DBI::errstr;
        $dbh->rollback();
        # stop processing IDs
        last Geneid;

    }

    #### Get polypeptide feature_id 
    my $pep_fids_ref = get_pep_fid(
                        $gene_id, $type_id_gene, $type_id_pseudo, $type_id_part, $type_id_derives,
                        $type_id_rrna, $type_id_trna
                        );

    #### Unfold multiple transcripts
    my @pep_fids = @{$pep_fids_ref};

    foreach my $pep_fid (@pep_fids){

        #### Unfold multiple PubMed IDs
        my @pubmed_ids = split(/\|/, $pubmed_ids{$gene_id});

        Pubmedid: foreach my $pubmed_id (@pubmed_ids){

            load_pubmed_ids_into_chado($gene_id, $pep_fid, $pubmed_id, $type_id_unfetched);
        }
    }
}

close $h_tag;
close $h_out;
close $h_drop;
close $h_temp;

#### Commit changes to the db
## uncomment when ready
$dbh->commit() unless(defined($errflag));

## close db connection
$dbh->disconnect();

################
sub get_pubmed_ids_from_tsv { 
    # stores pubmed_ids and gene_ids in a hash keys:ids, values:pubmed_ids 

    my $s_f_data = shift @_;
    my %s_pubmed_ids;

    open (my $h_data, "<", $s_f_data);

    # Get rid of the header line
    my $s_header = <$h_data>;

    while (<$h_data>){
        chomp;
        my $s_line = $_;
        my ($s_gene_id, $s_pubmed_id) = split (/\t/, $s_line);

        $s_pubmed_id =~ s/PMID://i;
        $s_pubmed_id =~ s/\s+//g;

        $s_pubmed_id = 'PMID:' . $s_pubmed_id;

        if(defined $s_pubmed_ids{$s_gene_id}){

            $s_pubmed_ids{$s_gene_id} =  $s_pubmed_ids{$s_gene_id} . '|' . $s_pubmed_id;

        }else{

            $s_pubmed_ids{$s_gene_id} = $s_pubmed_id;

        }
    }

    close $h_data;

    #### Return hash
    return \%s_pubmed_ids;

}

########

sub get_pep_fid {
    my ($s_gene_id, $s_type_id_gene, $s_type_id_pseudo, $s_type_id_part, $s_type_id_derives, @s_nc_types) =  @_;
    #my @s_nc_types = ($s_type_id_rrna, $s_type_id_trna ); # !00 For now only those two types in fungidb chado
    my ($s_gene_fid, $s_mRNA_fid);
    my @s_pep_fid;

    #### get gene feature id
    $s_gene_fid = $dbh->selectrow_array(
                                    $s_sql_get_gene_fid, undef, $s_gene_id, 
                                    $s_type_id_gene, $s_type_id_pseudo
                                    );

    if(defined $s_gene_fid){

        #### Get polypeptide feature id. Consider multiple transcripts 
        my @s_mRNA_fid = @{$dbh->selectcol_arrayref($s_sql_subject_id, undef, $s_gene_fid, $s_type_id_part)};

        foreach my $s_mRNA_fid (@s_mRNA_fid){
            # get polypeptide feature id for each transcript

            my $s_pep_fid=$dbh->selectrow_array($s_sql_subject_id, undef, $s_mRNA_fid, $s_type_id_derives);

            print "$s_gene_id\t$s_pep_fid\tPEP_FEATURE_ID\n";

            if(defined $s_pep_fid){ # check if there is a polypeptide (non coding genes)
                push(@s_pep_fid, $s_pep_fid);

            } else{ # check if the feature is part of a non coding gene

                my $s_nc_fid = $dbh->selectrow_array($s_sql_check_non_coding, undef, $s_mRNA_fid, @s_nc_types); 

                if (defined $s_nc_fid){

                    push(@s_pep_fid, $s_nc_fid); # add publication in the transcript feature

                }; # if not found will default to undefined @s_pep_fid and exit

            }
        }

    } else{

        print $h_drop "$s_gene_id\tID_NOT_FOUND\n";

    }


    unless( @s_pep_fid ){ # polypeptide feature_id not found
        #$errflag = 1; # Let the loading continue but log the missing gene/polypeptide ID
        print $h_drop "$s_gene_id\tPOLYPEPTIDE_NOT_FOUND\n";
        next Geneid; 
    }

    #### return polypeptide feature id
    return \@s_pep_fid;

}

########

sub load_pubmed_ids_into_chado { 
    ## Adds a PubMed ID to a polypeptide feature, checking that it does not exists

    my ($s_gene_id, $s_pep_fid, $s_pubmed_id, $s_type_id_unfetched) = @_; 

    ## Get pub_id (insert new one if necessary)
    my $s_pub_id = get_pub_id($s_pubmed_id, $s_type_id_unfetched);

    if(defined $s_pep_fid && defined $s_pub_id){

        ## Check if the pubmed reference is already added in the feature_pub table
        my $s_feature_pub_id = check_feature_pub($s_pep_fid, $s_pub_id);

        if(defined $s_feature_pub_id){ 
            # log and go to next reference
            
            print $h_drop "$s_gene_id\t$s_pubmed_id\tFPUB_ID_PRESENT\n";
            next Pubmedid; # Go to the next pubmed reference

        }else{

            ## Insert pubmed_id in the feature_pub table
            unless ($s_sql_add_feature_pub_id->execute($s_pep_fid, $s_pub_id) ){

                    # error inserting pub_id in the feature_pub table
                    $errflag = 1;
                    next Geneid;  # Rollback outside subroutine 
            }
        }

        print $h_out "$s_gene_id\t$s_pubmed_id\tADDED_FEATURE_PUB_ID\n";
        print $h_tag "$s_pep_fid\n";

    }
}

########
sub get_pub_id{
    # check if the PubMed ID entry is present in the pub table, insert it if it is not,
    # return the pub_id

    my ($s_pubmed_id, $s_type_id_unfetched) = @_;
    my $s_pub_id;
    my $s_obsolete = 'f';

    #### Check if the PubMed ID exists already in the pub table 

    $s_pub_id=$dbh->selectrow_array($s_sql_get_pub_id, undef, $s_pubmed_id, $s_type_id_unfetched);

    unless(defined $s_pub_id){

        unless($s_sql_add_pub_id->execute($s_pubmed_id, $s_type_id_unfetched, $s_obsolete)){

            ## error inserting
            print $h_drop "$s_pubmed_id\tERROR_INSERT_PUBMED_ID\n"; # log
            $errflag=1; # Rollback outside the subroutine
            next Geneid;
        }

        print $h_out "$s_pubmed_id\tADDED_PUBMED_ID\n";

        #### Get pub_id
        $s_pub_id = $dbh->selectrow_array($s_sql_get_pub_id, undef, $s_pubmed_id, $s_type_id_unfetched);
    }

    # return pub_id
    return $s_pub_id;
}

########
sub check_feature_pub{
    # Checks if a publication is already present in feature_pub table
    my ($s_pep_fid, $s_pub_id) = @_; 


    my $s_fpub_id=$dbh->selectrow_array($s_sql_get_fpub_id, undef, $s_pep_fid, $s_pub_id);

    return $s_fpub_id;
}
