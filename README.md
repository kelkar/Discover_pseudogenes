# Search_pseudogenes

## What this tool does:

(still working on this)

## Getting inputs ready

Before running this tool, you need to run exonerate, using protein sequences as queries and a genome as a target. Use the following exonerate parameters:

```
exonerate --percent 60 --model protein2genome --showquerygff yes --showtargetgff yes -q QUERY_PROTEINS.fasta -t TARGET_GENOME.fasta --ryo "RYO\t%qi\t%ti\t%ql\t%tl\t%qal\t%qab\t%qae\t%tal\t%tab\t%tae\t%et\t%ei\t%es\t%em\t%pi\t%ps\t%g\nTransitionStart\n%V{%Pqs\t%Pts\t%Pqb\t%Pqe\t%Ptb\t%Pte\t%Pn\t%Pl\n}TransitionEnd\nTargetSeq\n%qs\nAligned Sequences\n>Q\n%qas\n>T\n%tas\nCoding Sequences\n>Q\n%qcs\n>T\n%tcs\n" >  EXONERATE_OUTPUT.TXT
```

Then, process the output file to extract stop codon and frameshift information:

```
grep -P "Query:|Target:|(Query range:)|(Target range:)|(^[A-Z]\tTAA\t)|(\sframeshift\s)|(^[A-Z]\tTAG\t)|(^[A-Z]\tTGA\t)|(^\*)"  EXONERATE_OUTPUT.TXT   > EXONERATE_ERRORS.TXT
```

Also use this perl tool (https://github.com/EVidenceModeler/EVidenceModeler/blob/master/EvmUtils/misc/Exonerate_to_evm_gff3.pl) by EVM to convert exonerate output file to GFF3 format file.

```
perl Exonerate_to_evm_gff3.pl EXONERATE_OUTPUT.TXT > EXONERATE_OUTPUT.GFF3
```

## Usage

Finally, use all these input files to create a table of possible pseudogenization events:

```
perl tabulate_stops_frameshifts.pl QUERY_PROTEINS.fasta EXONERATE_OUTPUT.GFF3 EXONERATE_ERRORS.TXT > STOPS_FRAMESHIFTS.TXT
```
