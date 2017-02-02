#! /usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use DBI;
use IO::Prompt;

# Adds names to a list of Gene IDs
# Input file in tab format, with headers, first column gene IDs, second column gene names
# only one pair geneid:name per row
# Checks if the name is already present, also checks if the gene already has a name and
# in that case loads the new value as synonym (checks if is already a synonym too)

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
    -a|data_file <tab file with gene ID (first column) and name to load (second column)>
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

my $dbh = DBI->connect($dbi_connect, $dbuser, $dbpass,{RaiseError => 0, AutoCommit => 0}) or die "Can't connect to db $dbname!\n";

print "Success connecting to db $dbname!\n";

my $errflag;

#### Fungal instance
my $type_id_gene = 3596;
my $type_id_pseudo = 3229;
my $type_id_synonym = 1509;

#### Prepare SQL statements
## Check current feature name
my $s_sql_get_name = $dbh->prepare('SELECT name FROM feature WHERE feature_id=?');

## Add gene name to a feature
my $s_sql_add_name = $dbh->prepare('UPDATE feature SET name=? WHERE feature_id=?');

## Get feature id of gene or pseudogene by gene ID 
my $s_sql_get_gene_fid = $dbh->prepare('SELECT feature_id FROM feature WHERE uniquename=? AND (type_id=? OR type_id=?)');

## Get synonyms by gene ID
my $s_sql_get_feature_synonyms = $dbh->prepare('SELECT synonym.name FROM synonym join feature_synonym using (synonym_id) where feature_id=?');

## Get synonym_id 
my $s_sql_get_synonym_id = $dbh->prepare('SELECT synonym_id FROM synonym WHERE lower(name) = lower(?) AND type_id=?');

## Insert synonym in synonym collection
my $s_sql_add_synonym = $dbh->prepare('INSERT INTO synonym (name, type_id, synonym_sgml) VALUES (?,?,?)');

## Insert synonym in feature
my $s_sql_insert_feature_synonym = $dbh->prepare('INSERT INTO feature_synonym (synonym_id, feature_id, pub_id, is_current) VALUES (?,?,1,?)');

#### Get names into hashes key:gene_id, value:gene_name 
open (my $h_drop, ">>", $f_drop);
open (my $h_temp, ">>", $f_temp);

my $names_ref = get_names_from_tsv($f_data);
my %names = %{$names_ref};

## Check extraction

foreach my $id ( sort ( keys %names) ){

    print $h_temp "$id\t$names{$id}\n";

}

open (my $h_out, ">>", $f_out);
open (my $h_tag, ">>" ,$f_tag);

#### Load names into chado 
Geneid: foreach my $gene_id ( sort (keys %names)){
    
    if (defined $errflag){
        warn "Error adding name into feature: " . $DBI::errstr;
        $dbh->rollback();
        # stop processing IDs
        last Geneid;

    }

    #### Get gene feature_id
    my $gene_fid = get_gene_fid($gene_id,$type_id_gene,$type_id_pseudo);

    #### Unfold multiple names
    my @names = split(/\|/, $names{$gene_id});

    foreach my $name (@names){

        load_name_into_chado($gene_id, $gene_fid, $name, $type_id_synonym);
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
sub get_names_from_tsv { # stores names and gene_ids in a hash keys:ids, values:names 
    my $s_f_data = shift @_;
    my %s_names;

    open (my $h_data, "<", $s_f_data);

    # Get rid of the header line
    my $s_header = <$h_data>;

    while (<$h_data>){
        chomp;
        my $s_line = $_;
        my ($s_gene_id, $s_name) = split (/\t/, $s_line);

        if(defined $s_names{$s_gene_id}){

            $s_names{$s_gene_id} =  $s_names{$s_gene_id} . '|' . $s_name;

        }else{

            $s_names{$s_gene_id} = $s_name;

        }
    }

    close $h_data;

    #### Return hash
    return \%s_names;

}

sub load_name_into_chado {
    my ($s_gene_id, $s_gene_fid, $s_name, $s_type_id_synonym) = @_; 

    #### check if there is a name already
    my $s_curr_name = $dbh->selectrow_array($s_sql_get_name, undef, $s_gene_fid);

    ## if name exists, check if it is the same or other
    ## check if the name field has been initialized, but it is empty in the db
    ## (Artemis leaves the name field empty when a gene name is removed)

    if (defined $s_curr_name && $s_curr_name ne ''){

        if ( lc("$s_curr_name") eq lc("$s_name") ){ # if it is the same name just log
                
            print $h_drop "$s_gene_id\t$s_name\tNAME_PRESENT\n";

        }else{ # add as new name as synonym 

            add_synonym($s_gene_id, $s_gene_fid, $s_name, $s_type_id_synonym);

        }


    }else { # if not defined add new name

        unless( $s_sql_add_name->execute($s_name, $s_gene_fid) ){

                ## Error updating
                $errflag = 1;
                next Geneid; # Rollback outside subroutine
        }

        unless (defined $errflag){

            print $h_out "$s_gene_id\t$s_name\tADDED_NAME\n";
            print $h_tag "$s_gene_id\n";
        }
    }

}

####
sub get_gene_fid {
    my ($s_gene_id, $s_type_id_gene, $s_type_id_pseudo) = @_;

    #### get gene feature id
    my $s_gene_fid = $dbh->selectrow_array($s_sql_get_gene_fid, undef, $s_gene_id,
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

####
sub add_synonym { # adds a synonym to a gene_id checking if it already exists
    my ($s_gene_id, $s_gene_fid, $s_name, $s_type_id_synonym) = @_;
    my $s_present = 0;

    my $s_curr_synonyms_arrayref = $dbh->selectall_arrayref(
                                            $s_sql_get_feature_synonyms, undef, $s_gene_fid
                                            );

    # get number of rows (synonyms)
    my $s_nrows = scalar (@{$s_curr_synonyms_arrayref});

    unless ($s_nrows == 0){

        foreach my $s_curr_synonym_ref (@{$s_curr_synonyms_arrayref}){
            
            my $s_curr_synonym = @{$s_curr_synonym_ref}[0];

            # check if is the same synonym we want to add
            if ( lc($s_name) eq lc($s_curr_synonym)){

                print $h_drop "$s_gene_id\t$s_name\tFSYNONYM_PRESENT\n";
                $s_present = 1;

                last; # exits loop, fsynonym already present
            };

        }
    }

    unless ($s_present == 1) { # insert synonym unless it is already present

        my $s_synonym_id = get_synonym_id($s_gene_id, $s_gene_fid, $s_name,$s_type_id_synonym);
        my $s_current = 't';
        
        unless ($s_sql_insert_feature_synonym->execute($s_synonym_id, $s_gene_fid, $s_current )){

                # error inserting synonym
                $errflag = 1;
                next Geneid; 
        }

        print $h_out "$s_gene_id\t$s_name\tADDED_FSYNONYM\n";
        print $h_tag "$s_gene_fid\n";

    }
}

####
sub get_synonym_id {    
    # get the synonym_id needed for inserting into feature_synonym
    # create a new synonym entry if it does not exist in the db

    my ($s_gene_id, $s_gene_fid, $s_name, $s_type_id_synonym) = @_;

    my $s_synonym_id = $dbh->selectrow_array($s_sql_get_synonym_id, undef, 
                                            $s_name, $s_type_id_synonym
                                            );

    unless (defined $s_synonym_id) { 
    # add synonym to synonym db collection if not present already

        unless ($s_sql_add_synonym->execute($s_name, $s_type_id_synonym, $s_name)){ 
            
            # error inserting synonym?
            $errflag = 1;
            next Geneid; 
        }

        unless (defined $errflag){

            print $h_out "$s_gene_id\t$s_gene_fid\t$s_name\tADDED_SYNONYM\n";

        }

        $s_synonym_id = $dbh->selectrow_array($s_sql_get_synonym_id, undef,
                                            $s_name, $s_type_id_synonym
                                            );
    }

    # return synonym_id
    return $s_synonym_id;
}
