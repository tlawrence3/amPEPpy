# amPEPpy
## Installation
1. Download amPEPpy and training data using the below bash command or the zip link: 
```bash
git clone https://github.com/tlawrence3/amPEPpy.git
```
2. Change into the `amPEPpy` directory
```bash
cd amPEPpy
```
3. We recommend using [anaconda](https://www.anaconda.com/products/individual) and conda enviroments for installing and using amPEPpy. The following commands will create a conda environment with all the required packages and activate it.
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
Here are the steps to get amPEPpy up and running to classify protein sequences as AMP or nonAMP

### Training random forest classifier
Now that amPEPpy is installed we need to train the machine-learning algorithm. To do this we need a positive dataset (AMP sequences) and a negative dataset (nonAMP) sequences in fasta format. These are located in the `training_data` folder.

1. To train the machine-learning with the same settings as the manuscript using the below command.
```bash
ampep train -p training_data/M_model_train_AMP_sequence.numbered.fasta -n training_data/M_model_train_nonAMP_sequence.numbered.proplen.subsample.fasta --seed 2012
```

This should create a file named `amPEP.model` which contains the saved random forest classifier

### Classifying sequences using the trained classifier
Use the below command to classify amino acid sequences in fasta format using the trained random forest. As an example we will use our positive training data.

```bash
ampep predict -m amPEP.model -i training_data/M_model_train_AMP_sequence.numbered.fasta -o results.tsv --seed 2012
```

This should result in a file named `results.tsv` that contains the classification results for the positive dataset.
