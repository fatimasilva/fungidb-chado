#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use DBI;
use IO::Prompt;

# Adds previous IDs to a list of GeneIDs. 
# Input file in tab format, with headers, first column gene IDs, second column previous IDs
# only one pair geneid:previousid per row
# Checks if the previous ID is already present

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
    -a|data_file <tab file with gene ID (first column) and previousid to load (second column)>
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
my $f_tag = $name_out . "_gene_feature_ids_history.tsv";


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
my $type_id_previous = 2843;

#### Prepare SQL statements !00 needs checking

## Get feature id of gene or pseudogene by gene ID (uses type_id_gene and type_id_pseudo
my $s_sql_get_gene_fid = $dbh->prepare('SELECT feature_id FROM feature WHERE uniquename=? AND (type_id=? OR type_id=?)');

## Get previous ids by feature (gene) id (uses type_id_previous)
my $s_sql_get_previous_ids = $dbh->prepare('SELECT synonym.name FROM feature_synonym join synonym using (synonym_id) where feature_id=? and type_id=?');

## Get previous_fid, the identifier in the synonym table (uses type_id_previous)
my $s_sql_get_previous_fid = $dbh->prepare('SELECT synonym_id FROM synonym WHERE lower(name) = lower(?) AND type_id=?');

## Insert new previous id entry in the synonym table (uses type_id_previous)
my $s_sql_add_synonym = $dbh->prepare('INSERT INTO synonym (name, type_id, synonym_sgml) VALUES (?,?,?)')    ;

## Insert previous_fid in feature_synonym table
my $s_sql_add_previous_id = $dbh->prepare('INSERT INTO feature_synonym (synonym_id, feature_id, pub_id, is_current) VALUES (?,?,1,?)');

#### Get previous IDs into hashes key:gene_id, value:previous_id 
open (my $h_drop, ">>", $f_drop);
open (my $h_temp, ">>", $f_temp);

my $previous_ids_ref = get_previous_ids_from_tsv($f_data);
my %previous_ids = %{$previous_ids_ref};

## Check extraction

foreach my $id ( sort ( keys %previous_ids) ){

    print $h_temp "$id\t$previous_ids{$id}\n";

}

open (my $h_out, ">>", $f_out);
open (my $h_tag, ">>" ,$f_tag);

#### Load previous_ids into chado 
Geneid: foreach my $gene_id ( sort (keys %previous_ids)){
    
    if (defined $errflag){
        warn "Error adding previous ID into feature: " . $DBI::errstr;
        $dbh->rollback();
        # stop processing IDs
        last Geneid;

    }

    #### Get gene feature_id
    my $gene_fid = get_gene_fid($gene_id, $type_id_gene, $type_id_pseudo);

    #### Unfold multiple names
    my @previous_ids = split(/\|/, $previous_ids{$gene_id});

    foreach my $previous_id (@previous_ids){

        load_previous_ids_into_chado($gene_id, $gene_fid, $previous_id, $type_id_previous);
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
sub get_previous_ids_from_tsv { 
    # stores previous_ids and gene_ids in a hash keys:ids, values:previous_ids 

    my $s_f_data = shift @_;
    my %s_previous_ids;

    open (my $h_data, "<", $s_f_data);

    # Get rid of the header line
    my $s_header = <$h_data>;

    while (<$h_data>){
        chomp;
        my $s_line = $_;
        my ($s_gene_id, $s_previous_id) = split (/\t/, $s_line);

        if(defined $s_previous_ids{$s_gene_id}){

            $s_previous_ids{$s_gene_id} =  $s_previous_ids{$s_gene_id} . '|' . $s_previous_id;

        }else{

            $s_previous_ids{$s_gene_id} = $s_previous_id;

        }
    }

    close $h_data;

    #### Return hash
    return \%s_previous_ids;

}

########

sub get_gene_fid {
    my ($s_gene_id, $s_type_id_gene, $s_type_id_pseudo) =  @_;

    #### get gene feature id
    my $s_gene_fid = $dbh->selectrow_array(
                                    $s_sql_get_gene_fid, undef, $s_gene_id, 
                                    $s_type_id_gene, $s_type_id_pseudo
                                    );

    unless(defined $s_gene_fid){ # gene id not found
        #$errflag = 1; # Let the loading continue but log the missing geneID
        print $h_drop "$s_gene_id\tNOT_FOUND\n";
        next Geneid; 
    }

    #### return gene_fid
    return $s_gene_fid;

}

########

sub load_previous_ids_into_chado { 
    ## Adds a previous_systematic_id to a gene feature, checking that it does not exists

    my ($s_gene_id, $s_gene_fid, $s_previous_id, $s_type_id_previous) = @_; 

    my $s_present = 0;
    
    ## Get a list of previous ids present in that feature
    my $s_curr_previous_arrayref = $dbh->selectall_arrayref(
                                            $s_sql_get_previous_ids, undef, $s_gene_fid, 
                                            $s_type_id_previous)
                                            ;

    # get number of rows (previous ids)
    my $s_nrows = scalar (@{$s_curr_previous_arrayref});

    unless ($s_nrows == 0){ # skip if there are no previous ids

        foreach my $s_curr_previous_ref (@{$s_curr_previous_arrayref}){ # unfold sql results
            
            my $s_curr_previous = @{$s_curr_previous_ref}[0];

            # check if is the same previous_id we want to add
            if ( lc($s_previous_id) eq lc($s_curr_previous)){

                print $h_drop "$s_gene_id\t$s_previous_id\tPREVIOUS_ID_PRESENT\n";
                $s_present = 1;

                last; # exits loop, previous_id already present
            };

        }
    }

    unless ($s_present == 1){ # insert previous id unless it is already present

        ## Get the identifier of the previous id to be inserted (synonym_id in the db) 
        ## if it does not exist in the db create an entry in the synonym table

        my $s_previous_fid = get_previous_fid(
                                $s_gene_id, $s_gene_fid, $s_previous_id, $s_type_id_previous);
        my $s_current = 't';
        
        ## Insert previous_fid in the feature_synonym table
        unless ($s_sql_add_previous_id->execute($s_previous_fid, $s_gene_fid, $s_current) ){

                # error inserting previous_fid
                $errflag = 1;
                next Geneid;  # Rollback outside subroutine 
        }

        print $h_out "$s_gene_id\t$s_previous_id\tADDED_PREVIOUS_ID\n";
        print $h_tag "$s_gene_fid\n";

    }
}

########

sub get_previous_fid {  
    # gets the previous_fid (db internal number) needed for inserting into feature_synonym
    # creates a new previous_id entry in the synonym table if it does not exist

    my ($s_gene_id, $s_gene_fid, $s_previous_id, $s_type_id_previous) = @_;

    ## Check if the previous_id has an entry in the database

    my $s_previous_fid = $dbh->selectrow_array(
                                $s_sql_get_previous_fid, undef, $s_previous_id, $s_type_id_previous
                                );

    unless (defined $s_previous_fid) { 
        # add previous_id to synonym db collection if not present already

        unless ($s_sql_add_synonym->execute($s_previous_id, $s_type_id_previous, $s_previous_id)){ 
            
            # error inserting synonym?
            $errflag = 1;
            next Geneid; 
        }

        unless (defined $errflag){

            print $h_out "$s_gene_id\t$s_gene_fid\t$s_previous_id\tADDED_PREVIOUS_ID_ENTRY\n";

        }

        $s_previous_fid = $dbh->selectrow_array(
                            $s_sql_get_previous_fid, undef, $s_previous_id, $s_type_id_previous);

    }

    # return previous_fid
    return $s_previous_fid;
}
