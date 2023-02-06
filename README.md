# General purpose

This pipeline can be used as template for novel pipelines that are build in the future. The structure of this README is also very silimar to the structure of the other pipelines so feel free to modify this accordingly. 
---

## How to use

If this is the first time that you are running this pipeline, please look at the read me in the [home page](https://github.com/AG-Walz#how-to-use---first-timers-edition). This would save time from both parties (yours and the bioinformaticians).

<details markdown="1">
<summary>Input setup</summary>

### Input file

Here, you can write the purpose of the input file and what it needs to contain. Furthermore, it is of interest to provide an example.


---

### File structure

When everything is set up properly, then your input directory would look like the following tree.

```bash
inputDir
  |_ data
  |_ data
```


---

</details>

<details markdown="1">
<summary>Output</summary>

### Output folder structure

When everything goes well, a series of files is created in the defined results folder.
During the rest of this section, we go over the different folders and the output that is reported in them.

```bash
outputDir
  |_ content
```

---

</details>

<details markdown="1">
<summary>Frequently reported problems</summary>

### Libraries

The most important thing is to check if all of the different libraries are there. When you stumble accross the following piece of code pay close attention!

```r
################################################################################
###                             Load libraries                               ###
################################################################################
## Load the libraries
necessaryLibs <- c("example")

## Load the necessary libraries()
invisible(lapply(necessaryLibs, library, character.only = T))
```

All of the libraries in the list of `necessaryLibs` needs to be installed since we need these packages for some of the functionalities to run. It could happen that you did not install a package yet, then you will get the following error:

```r
library(example)
Error in library(example) : there is no package called ‘example’
```

To resolve this problem we have to install the packages, which we do using the following command `install.packages("example")`.
If that does not work, then we can use `BiocManager::install("example")`.
Make sure that you put the qoutation marks around the package of interest.

If everything is installed, and `invisible(lapply(necessaryLibs, library, character.only = T))` runs without errors, then we can finally run the rest of the code

> Note, it could be that you see red text, and **error** is only an error when it is mentioned explicity!

---

</details>
