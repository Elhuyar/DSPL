#!/usr/bin/perl

my $usage = <<"_USAGE_"; 

    separateCorpubyPol.pl

    Functionality: this script divides a corpus according to the polarity of its documents. Documents are limited by <doc>text</doc> tags. Polarity is deined by the attribute <doc pol="(pos|neg|neu)">. 

    Parameters:

        1.- Corpus file name.

    Call Example: % perl separateCorpusbyPol /path/to/corpus/corpus.txt

    Output: One file per each of the polarities found in the input file, containing the documents with that polarity. For example:

               /path/to/corpus/corpus_pos.txt
               /path/to/corpus/corpus_neg.txt
               /path/to/corpus/corpus_neu.txt
                        ...  

    To see documentation and usage, co:

    separateCorpubyPol.pl

_USAGE_

my $file=shift;

($basename, @exts) = split /\./, $file;
my $ext=join('.',@exts);

open(INFILE, "<$file") || die "$usage\n ERROR: The corpus file was not specified or could not be opened.\n";

while ($l=<INFILE>)
{
    #<doc> irteerako fitxategia aldatu kategoriaren arabera
    if ($l =~ /<doc .*pol="([^\"]*)"/)
    {
	$pol=$1;	
	open(OUTFILE, ">>$basename\_$pol.$ext");
	print OUTFILE "$l";
	
    }
    elsif ($l =~ /<\/doc>/)
    {
	print OUTFILE "</doc>\n";
	close(OUTFILE);
    }
    else
    {
	print OUTFILE "$l";
    }
    
}
close(OUTFILE);
close(INFILE);
