# Description

​

## Introduction

The project idea is to identify recombination events among herpes viruses. Recombination is a primary source of genetic variation among certain types of viruses, such as Herpesviridae. Recombination occurs when two species of viruses infect the same cell. With two viral genomes in the same cell, there is an opportunity for the genomes to exchange and join regions of their genetic material. The exact details of how recombination events occur are complex, and the identification of recombination events is not straightforward. The phenomenon is known since 1955 and is believed to be historically important for the evolution of herpesviruses, both within strains and between species.

## Assumptions and intuition

I assume that the nucleotide composition in regions where a foreign stretch of DNA has been inserted is divergent from the rest of the genome. This makes intuitive sense since it is established that viruses have specific k-mer compositions. Therefore, it should be possible to identify recombination events by constructing a statistical model of the genome and scanning for parts of the genome where there is a large degree of dissimilarity to the model. However, such regions can also have other explanations, e.g. if the regions are conserved within entire families.

## Results

I've put all of the plots in a [notebook](Recombination.ipynb), together with visualisation code. Most of the heavy lifting is done in the [util](util) module.

My first idea was to train an HMM on the full genome sequence and then compute the log-likelihood of subsequences of the genome under the HMM. However, this turns out to give the same result as computing the GC-content of the regions. As such, the HMM is not able to identify any specific regions within the genome to a greater extent than just the GC content would.

I then tried to compute the 4-mer (among others) frequencies (normalised by the nucleotide frequencies) and compared these to the subsequences of the genomes. The intuition is that we would identify subsequences that are more similar to a reference genome. Plotting these results, however, does not indicate that we can use this for identification of recombination events. A statistic test here identifies regions where two genomes are statistically different, while this is the opposite test from what we are after it illustrates that not many regions are statistically different. Thus, this method is not sensitive enough.

Finally, I turned to raw long k-mers. The intuition here is that in the event of a recombination event, the two genomes should share long k-mers. However, long k-mers can also be expected to exist in viruses with a single common ancestor. On the other hand, if two viral species share a long k-mer while that k-mer is not shared with a putative ancestor (e.g. if it is not conserved within the genus/subfamily), we could hypothesise that the k-mer is a result of a recombination event.

## Conclusion

Recombination events are difficult to identify using probabilistic and alignment-free models. The most promising idea tested here consisted of identifying shared long k-mers. However, the origin of the shared k-mers needs further analysis. They could just as well be a result of random point-mutations or shared ancestors rather than recombination.

### Other approaches

We also discussed using profile HMMer. These models may have been more successful. However, I felt that by starting with a multiple sequence alignment (which is the first step of the HMMer profile HMMs) would by definition require that the recombination events occurred at the same position throughout the viruses, which of course is not the case. However, if we were searching for a specific recombination event (previously identified) this would probably have worked well.

Sources:

[Recombination in viruses: Mechanisms, methods of study, and evolutionary consequences](https://www.sciencedirect.com/science/article/pii/S156713481400478X)

[Divergence and genotyping of human α-herpesviruses: An overview](https://www.sciencedirect.com/science/article/pii/S1567134809001993)
