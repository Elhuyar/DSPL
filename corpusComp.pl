#!/usr/bin/perl


# command to get a frequnecy list out of a freeling tagged file
## grep -v "^</\?doc" corpus.fling.tagged | grep -v "^\s*$" | cut -f2,3 -d' ' | sort | uniq -c | sed -e "s/^\s\+//" | sort -nr > corpus.freq
# command to get a frequnecy list out of a freeling tagged file (including sense information) 
# grep -v "^</\?doc" corpus.fling.tagged | grep -v "^\s*$" | cut -f2,5 -d' ' | sed -e "s/\(-[arvn]\):.*$/\1/" | sort | uniq -c | sed -e "s/^\s\+//" | sort -nr > corpus.freq

my $usage = <<"_USAGE_"; 
    
    corpusComp.pl  

    Functionality: Compares two corpora, extracting the most salient words of the first corpus with respect to the other, according to different statistical measures. This is useful for example to find domain terminology by comparing a domain corpus vs. a general domain corpus, or finding subjective words by comparing a corpus of subjective documents vs. and objective corpus. Another example is to extract polarity words by comparing positive and negative review corpora. 

    Frequency Lists: The scripts takes as input the two lists of lemma frequencies corresponding to each corpus. It is up to the user to provide those lists, but the following commandline sequence creates such list for a corpus tagged with FreeLing or TreeTagger:

    % grep -v "^</\?doc" corpus.fling.tagged | grep -v "^\\s*\$" | cut -f2,3 -d\' \' | sort | uniq -c | sed -e "s/^\\s\\+//" | sort -nr > corpus.freq

    Parameters:
        1. domain corpus (-c | --corpus= ): frequency list of the lemmas of the domain / specific corpus we want to extract keywords from. Format is tabulated format "frequency \t lemma".    
        2. reference corpus  (-r | --reference= ): frequency list of the lemmas of the reference corpus used in the comparison. Format is tabulated format "frequency \t lemma". 
        3. association measure(-m | --measure= ): association measure to be used on comparing corpora. Default value is llr. possible values: (ll|log-odds|odds|pmi|fisher|diff).
        4. -v | --verbose : print extra info.
       
    Output: 
          Standard output.  

    Call Examples:
           $ perl corpuscomp -c inputCorpus1.txt -r inputCorpus2.txt  ---- results to standard output
           $ perl corpuscomp -c inputCorpus1.txt -r inputCorpus2.txtCorpus.txt -m pmi > outputFile.txt  ---- results to file, association measure pmi (point mutual information).


    To see documentation and usage, type:

    corpusComp.pl --help


_USAGE_

#############################################################
#
#          Main Program
#
#############################################################

use Getopt::Long;
use Switch;
require Text::NSP::Measures::2D::MI::ll;
require Text::NSP::Measures::2D::Fisher::left;
require Text::NSP::Measures::2D::MI::pmi;
require Text::NSP::Measures::2D::odds;


#parameter initialization
my $synset=0;
my $verbose=0;
my $refFile="";
my $corpFile="";
my $measure="ll";
my %params;
(GetOptions (\%params,
	    "corpus|c=s"=>\$corpFile,
	    "reference|r=s"=>\$refFile,
	    "measure|m=s"=>\$measure,
	    "verbose|v", 
	    "help|h")) || die "ERROR: Illegal arguments or parameters: @ARGV\n" unless ($#ARGV < 0);

# parameter reading
if (defined $params{"help"} || defined $params{"h"})
{
    die "$usage\n";
}

if ($corpFile eq "" || $refFile eq "")
{
    die "Either specific or reference corpora are missing, both are required to run this script. \n $usage \n";
}

#word vs. synset based lexicon.
if ((! defined $params{"measure"}) && (! defined $params{"m"}))
{
    print STDERR "No measure provided Log Likelyhood Ratio (LLR) will be used\n";
}

#print extra info.
if ($params{"verbose"} || $params{"v"})
{
    $verbose=1;
}

# measure case insensitive.
$measure = lc($measure);


# global variables
my ($sum_corp_wc, %measure_values, %ExpectedFreqs, $hitzaFreq_ref, $hitzaFreq_spec, $measure_score);


# load both corpora
# For the moment size of the corpora is computed by suming up the frequencies of the lemmas. This is right if the frequency list includes all tokens (including function words and other symbols).
# IMPORTANT:  only nouns (N), adjective (JJ), verbs (V) and adverbs (RB) are taken as keywords, (Treebank tagset is used).
(my $specCorpus, my  $specCorpus_wc)=&load_freqList($corpFile);
$kk= scalar (keys (%$specCorpus));
print STDERR "corpusComp.pl: specific corpus loaded, size: $specCorpus_wc ; contentwords: $kk\n";
(my $refCorpus, my $refCorpus_wc)=&load_freqList($refFile);
$kk= scalar (keys (%$refCorpus));
print STDERR "corpusComp.pl: reference corpus loaded, size: $refCorpus_wc ; contentwords: $kk\n";


#Bi corpusetako hitz kopurua gorde
$sum_corp_wc = $specCorpus_wc + $refCorpus_wc; 
#print "Both corposa have $sum_corp_wc words altogether.\n";

for my $hitza (keys %$specCorpus)
{
    if(defined $refCorpus->{$hitza})
    {	
	$hitzaFreq_ref = $refCorpus->{$hitza};
    }
    else
    {
	$hitzaFreq_ref = 0;
    }

    $hitzaFreq_spec = $specCorpus->{$hitza};
    
    my $sum_hitzaFreqs = $hitzaFreq_spec + $hitzaFreq_ref;    

    #my$E1 = $specCorpus_wc *($sum_hitzaFreqs)/$sum_corp_wc;
    #my$E2 = $refCorpus_wc *($sum_hitzaFreqs)/$sum_corp_wc;

################################################
# 
#             |corpus1 | corpus2  | total     |
#  __________ |_______ |_________ |__________ |
#    w1       |  a=n11 |  b=n12   |  a+b=n1p  |
#  __________ |_______ |_________ |__________ |
#    not w1   |  c=n21 |  d=n22   |  c+d=n2p  |
#  __________ |_______ |_________ |__________ |
#             | a+c=np1| b+d=np2  |a+b+c+d=npp|
#             |        |          |           |


#LLR balioaren kalkulua

    my $n11=$hitzaFreq_spec;
    my $n1p=$sum_hitzaFreqs;
    my $np1=$specCorpus_wc;
    my $npp=$sum_corp_wc;

    
    my $E11 = ($n1p*$np1)/$npp;
   #if ($n11 < 2)
   #{	
	#next;
   #}
   if($measure eq "ll")
    { 
	$measure_score = Text::NSP::Measures::2D::MI::ll::calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);

	if($n11 < $E11){

	   $measure_score = -$measure_score;
	} 
    }
    elsif($measure eq "fisher")
    {
	$measure_score = Text::NSP::Measures::2D::Fisher::left::calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
    }  
    elsif($measure eq "pmi")
    {	
	$measure_score = Text::NSP::Measures::2D::MI::pmi::calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
    }
    elsif($measure eq "odds")
    {   
        $measure_score = Text::NSP::Measures::2D::odds::calculateStatistic(n11=>$n11,n1p=>$n1p,np1=>$np1,npp=>$npp);
    }
    elsif($measure eq "log-odds")
    {   
        #$n11 = $n11 + 0.5; 
        #$n1p = $n1p + 0.5;
        #$np1 = $np1 + 0.5; 
        #$npp = $npp + 0.5; 
	
	$n11 = $n11 + 0.5; 
        $n1p = $hitzaFreq_ref + 0.5;
        $np1 = $specCorpus_wc-$hitzaFreq_spec + 0.5; 
        $npp = $refCorpus_wc-$hitzaFreq_ref + 0.5;          
        
	my$goiko_zatia = $n11 * $npp; 
        
        my$beheko_zatia = $n1p * $np1; 
        
        my$emaitza = $goiko_zatia / $beheko_zatia; 

	# log = ln
        $measure_score = log($emaitza);
	# log base 10
        #$measure_score = log($emaitza)/log(10);	
    }
    elsif($measure eq "diff")
    {
	
	#warn "$hitza: $n11 -$np1 - $refCorpus_wc \n";

	# Maiztasun normalizatuen kalkulua (hitzaren maiztasuna 1M hitzetan)
	$nf1 = (($n11+1)*1000000)/$np1;
	$nf2 = (($hitzaFreq_ref+1)*1000000)/$refCorpus_wc;
	
        # Maiztasun normalizatuen ordez erlatiboak erabili nahi izanez gero:
        #$nf1 = $n11/$np1;
	#$nf2 = $hitzaFreq_ref/$refCorpus_wc;

	if ($nf2 == 0)
	{
	    $nf2 = 0.0001
	}
	
	#if (($nf1 > $nf2) && ($nf1 > 2) && ($nf2 > 2))
	# Maks & Voessen (2012)-ren arabera, min_freq(w)>2 eta freq=absolute frequency
	#if (($nf1 > $nf2)) #&& ($n11 > 2) && ($hitzaFreq_ref > 2))
	#if (($nf1 > $nf2)) #&& ($n11 > 2) && ($hitzaFreq_ref > 2))
	#{
	    $measure_score=(($nf1-$nf2)*100)/$nf2;
	#}
	#else
	#{
	#    $measure_score="eee"
	#}
    }
    if ($measure_score ne "eee")
    {
	$measure_values{$hitza} = $measure_score;   
	$ExpectedFreqs{$hitza} = $E11;   
    }
}


foreach $value (sort {$measure_values{$b} <=> $measure_values{$a} }
           keys %measure_values)
{
    if (! defined $refCorpus->{$value})
    {
	$refCorpus->{$value}=0;
    }
    if ($measure ne "diff")
    {
	print "$value\t$measure_values{$value}\t$specCorpus->{$value}\t$refCorpus->{$value}\t$ExpectedFreqs{$value}\n";  
    }   
    else
    {
	print "$value\t$measure_values{$value}\t$specCorpus->{$value}\t$refCorpus->{$value}\n";     
    }
}



#############################################################
#
#          Functions 
#
#############################################################


sub load_freqList()
{
    my $path=shift;

    my %list=();
    my $size=0;
    open(LEX, "<$path") || die ("corpusComp.pl: Could not load the polarity lexicon \'$path\'. It does not exist or it can't be opened. \n");
    while ($l=<LEX>)
    {
	chomp $l;
	# comment line
	#if (($l =~ /^# [A-Z]/) || ($l =~ /^\s*$/))
	#{
	#    next;
	#}
	(my $freq, my $lemma, my $pos, my @sobra) = split /\s+/, $l; 


	if ($pos=~/^(JJ|RB|V|NN)$/) # ($pos=~/^(Varn]$/)|| ($pos=~/^(Varn]$/))
	{
	    # normalize pos to single letter wordnet notation
	    $pos=~s/^V.*/v/;
	    $pos=~s/^RB/r/;
	    $pos=~s/^JJ/a/;
	    $pos=~s/^N.*/n/;

	    $list{"$lemma\_$pos"}=$freq;
	}
	$size+=$freq;
    }
    close(LEX);
    
    return (\%list, $size);

}
