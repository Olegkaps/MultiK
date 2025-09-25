MultiK is a tool for handling multiply mapped reads in DNA-DNA / DNA-RNA interactome analysis.

It is a powered version of mHi-C, providing the ability to use RNA annotation as bins and account for contact density.

Also, it works faster than mHi-C.

To use the default version, use Docker (you must have input files written in pipeline.sh).

```

docker build . -t multik

docker run multik

```

As written in Dockerfile, it runs ./compile.sh and ./pipeline.sh.

./DevideUni is needed in DNA-RNA mode, as well as genomic annotation and mapping of chromosome IDs in annotation to names in the input file.

In DNA-DNA mode skip ./DevideUni and ./MergePriorProbs steps; just use ./PriorC to parsed uni contacts.

Option -d of ./PriorC creates a file of DNA contact density, which is used in ./Multi2Uni.

Extra information is provided in the helpers of programs:

```

<programm name> --help

```

Thanks to the Mironov lab, this project is frozen and will never be updated. It has defects built in by design and had no real review at any stage.

Due to poor code quality and system design standards, every added feature made this project more useless and hard to support.
This project started from fake data showing that mHi-C works for an extremely long time. And from when falsification is visible, the program has no future.
