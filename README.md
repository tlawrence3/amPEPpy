# amPEPpy: An antimicrobial peptide prediction tool

## About
amPEPpy is a Python 3 application that implements a random forest classifier for predicting antimicrobial peptide sequences. amPEPpy has improved portability, increased accuracy relative to similar methods, and includes utilities for easily training and optimizing random forest classifiers on novel training data.
## Table of Contents


 * [Install](#install)
 * [Quickstart tutorial](#quickstart-tutorial)
     * [Training random forest classifier](#training-random-forest-classifier)
     * [Classifying sequences using the trained classifier](#classifying-sequences-using-the-trained-classifier)
 * [Optimizing the random forest classifier, calculating feature importance, and feature selection](#optimizing-the-random-forest-classifier-calculating-feature-importance-and-feature-selection)
     * [Optimizing the number of decision trees within the random forest classifier](#optimizing-the-number-of-decision-trees-within-the-random-forest-classifier)
     * [Feature importance](#feature-importance)
     * [Excluding features from training and classifiying](#excluding-features-from-training-and-classifiying)
 * [Citing](#citing)
 * [Appendix](#appendix)
     * [Global options](#global-options)
     * [Train](#train)
     * [Predict](#predict)


## Install
1. Download amPEPpy and training data using the below bash command or the zip link: 
```bash
git clone https://github.com/tlawrence3/amPEPpy.git
```
2. Change into the `amPEPpy` directory
```bash
cd amPEPpy
```
3. We recommend using [anaconda](https://www.anaconda.com/products/individual) and conda enviroments for installing and using amPEPpy. For MacOSX and Linux users we recommend installing the command line tools for anaconda. Windows users don't have this option and will need to open a powershell from the anaconda GUI to install and use `amPEPpy`. The following commands will create a conda environment with all the required packages and activate it.
```bash
conda create -n amPEP python=3.8 pandas numpy biopython scikit-learn
conda activate amPEP
```
4. Now we can install amPEPpy and test the installation with the below commands:
```bash
python setup.py install
ampep -h
```

## Quickstart tutorial
Here are minimal steps required to get `amPEPpy` up and running to classify protein sequences as AMP or nonAMP

### Training random forest classifier
Now that amPEPpy is installed we need to train the machine-learning algorithm. To do this we need a positive dataset (AMP sequences) and a negative dataset (nonAMP) sequences in fasta format. These are located in the `training_data` folder.

To train the random forest classifier with the same settings as the manuscript use the below command.
```bash
ampep train -p training_data/M_model_train_AMP_sequence.numbered.fasta -n training_data/M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta --seed 2012
```

This should create a file named `amPEP.model` which contains the saved random forest classifier. The `-p` flag is the path to your fasta file containing AMP sequences, `-n` flag is the path to a fasta file containing nonAMP sequences, and `--seed` may be provided for reproducibility and was set to `2012` for all of the analyses for the manuscript.

### Classifying sequences using the trained classifier
Use the below command to classify amino acid sequences in fasta format using the trained random forest. As an example we will use our positive training data.

```bash
ampep predict -m amPEP.model -i training_data/M_model_train_AMP_sequence.numbered.fasta -o results.tsv --seed 2012
```

This should result in a file named `results.tsv` that contains the classification results for the positive dataset. The `-m` flag is the path to the model file created during training, `-i` flag is the path to a fasta file containing amino acid sequences to classifier, and the `-o` flag. 

## Tutorial for optimizing the random forest classifier on novel training data
After finishing the quickstart tutorial the steps belows will take you through the steps we used to optimize our random forest classifier. These steps can be used to optimize a new classifier on novel training data. 

### Optimizing the number of decision trees within the random forest classifier
First we need to determine the number of decision trees produces the lowest [out-of-bag error](https://towardsdatascience.com/what-is-out-of-bag-oob-score-in-random-forest-a7fa23d710). The below command will calculate out-of-bag error for random forest classifiers containing between 23 and 175 decision trees saving the results to `oob_error_results.tsv`.
```bash
ampep train --test-trees --min-trees 23 --max-trees 175 \ 
-p training_data/M_model_train_AMP_sequence.numbered.fasta \
-n training_data/M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta --seed 2012 > oob_error_results.tsv
```

### Calculating feature importance
### Excluding features from training and classifiying

## Citing
## Appendix
### Global options
`-v --version`
Show program's version number and exit

`-h --help`
Print help message and exit

`--seed`
Seed for random processes. This allows reproducibility of random forest training. Default is a random seed number.

`-t --num-processes`
Number of processor cores to use for training and classification. Default is the number of system cores reported by the OS

`-d --drop-features`
Text file containing list of features to drop during random forest training and classification. File must have one feature per line

### Train
`-p --positive`
Fasta file containing amino acid sequences of known antimicrobial peptides.

`-n --negative`
Fasta file containing amino acid sequences of non-antimicrobial peptides.

`--test-trees`
Test accuracy of random forest classifier using different numbers of classifiers evaluted using out of bag error. `--min-trees` and `--max-trees` control the range of classifiers tested. Results are printed to stdout.

`--min-trees`
Minimum number of classifiers within random forest classifier when evaluating out of bag error. Default is 23.

`--max-trees`
Minimum number of classifiers within random forest classifier when evaluating out of bag error. Default is 175.

`--num-trees`
Number of classifers used to train random forest classifier. Default is 160 which was shown to produce the lowest out of bag error on training data.

`--feature-importance`
Test feature importance using the drop feature method. Results are written to a csv file in the current directory.

### Predict
`-i --input-sequences`
Fasta file containing amino acid sequences to be classified.

`-o --output-file`
Output file to save classification results. Default is stdout.

`-m --model`
Random forest model to use for classification. Random forest model is trained and outputted using the train module. Default is amPEP.model in the current directory.
