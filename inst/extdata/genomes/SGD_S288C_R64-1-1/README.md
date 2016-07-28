# Reference Genome SGD_S288C_R64-1-1

This directory contains information about the SGD reference genome SGD_S288C_R64-1-1.

Reference sequence information is contained in the following file:

## seqtab.csv

This CSV file contains information about reference sequences under the headings:

- `seqids`: normalised reference sequence identifiers (e.g. `04`).
- `seqnames`: canonical sequence names (e.g. `IV`).
- `seqlengths`: reference sequence physical map lengths (in base pairs). These
                are taken directly from the reference FASTA file of the
                [SGD_S288C_R64-1-1 assembly](http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz) of
                SGD (Cherry et al. 2012).
- `maplengths`: reference sequence genetic map lengths (in centiMorgans). Map
                lengths for the nuclear chromosomes are calculated from the
                per-chromosome genetic/physical distance ratios for the combined
                physical and genetic maps of *Saccharomyces cerevisiae* as given in the 
                [SGD Wiki](http://wiki.yeastgenome.org/index.php/Combined_Physical_and_Genetic_Maps_of_S._cerevisiae).
                The genetic map length of the mitochondrial chromosome ('17') is
                calculated from the genome-wide mean recombination rate of
                0.333 cM/kb (Petes et al. 1991, cited in Petes 2001).
- `isCircular`: `TRUE` if reference sequence is circular; `FALSE` otherwise.
- `genome`: contains the name of the reference genome assembly of which the
            given sequence is a part.

## References

Cherry JM, Hong EL, Amundsen C, Balakrishnan R, Binkley G, Chan ET, Christie
KR, Costanzo MC, Dwight SS, Engel SR, Fisk DG, Hirschman JE, Hitz BC, Karra K,
Krieger CJ, Miyasato SR, Nash RS, Park J, Skrzypek MS, Simison M, Weng S, Wong
ED. (2012) Saccharomyces Genome Database: the genomics resource of budding yeast.
Nucleic Acids Research. 40(Database issue):D700-5.

Petes TD, Malone RE, Symington LS (1991) 'Recombination in yeast.' in The
molecular and cellular biology of the yeast Saccharomyces. [Vol.1], Genome
dynamics, protein synthesis and energetics. Cold Spring Harbor : Cold Spring
Harbor Laboratory Press.

Petes TD. (2001) Meiotic recombination hot spots and cold spots.
Nature Reviews Genetics. 2(5):360-9.

