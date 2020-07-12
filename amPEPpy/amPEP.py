import sys
import os
import argparse
import random
import pickle
import numpy as np
import pandas as pd
import sklearn.utils
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import plot_roc_curve
from sklearn import metrics
from Bio import SeqIO
from amPEPpy._version import __version__

def main():
    parser = argparse.ArgumentParser(prog="ampep",
                                     description="amPEPpy is a python implementation of the amPEP random forest approach for identifying antimicrobial peptides. In addition to the original implementation, we have included a commandline interface, improved performance to scale to genomic data, and improved feature selection procedures",
                                     epilog="Please cite Lawrence et al. (2020) AmPEPpy: Antimicrobial Peptide prediction in python and Bhadra et al. (2018) AmPEP: Sequence-based prediction of antimicrobial peptides using distribution patterns of amino acid properties and random forest.",
                                     add_help=False)
    parser.add_argument('-v', '--version', action='version', version="v{}".format(__version__))
    parser.add_argument("--seed", type=int, default=random.SystemRandom().randint(0, 1000000), dest="seed",
                        help="Seed for random processes. This allows reproducibility of random forest training. Default is a random seed number.")
    parser.add_argument("-t", "--num-processes", type=int, default=os.cpu_count(), dest="num_processes",
                        help="Number of processor cores to use for training and classification. Default is the number of system cores reported by the OS")
    parser.add_argument("-d", "--drop-features", dest="drop_feature", help="Text file containing list of features to drop during random forest training and classification. File must have one feature per line")
    subparsers = parser.add_subparsers(help="sub-command help")

    parser_train = subparsers.add_parser("train", parents=[parser], help="Module to train random forest classifier to identify antimicrobial peptides from amino acid sequence data.")
    requiredNamed_train = parser_train.add_argument_group('required arguments')
    requiredNamed_train.add_argument("-p", "--positive", help="Fasta file containing amino acid sequences of known antimicrobial peptides.", required=True)
    requiredNamed_train.add_argument("-n", "--negative", help="Fasta file containing amino acid sequences of non-antimicrobial peptides.", required=True)
    parser_train.add_argument("--test-trees", action="store_true", dest="tree_test", help="Test accuracy of random forest classifier using different numbers of classifiers evaluted using  out of bag error. --min-trees and --max-trees control the range of classifiers tested. Results are printed to stdout.")
    parser_train.add_argument("--min-trees", type=int, default=23, dest="min_tree", help="Minimum number of classifiers within random forest classifier when evaluating out of bag error. Default is 23.")
    parser_train.add_argument("--max-trees", type=int, default=175, dest="max_tree", help="Minimum number of classifiers within random forest classifier when evaluating out of bag error. Default is 175")
    parser_train.add_argument("--num-trees", type=int, default=160, dest="num_trees", help="Number of classifers used to train random forest classifier. Default is 160 which was shown to produce the lowest out of bag error on training data.")
    parser_train.add_argument("--feature-importance", action="store_true", dest="feature_importance", help="Test feature importance using the drop feature method. Results are written to a csv file in the current directory.")
    parser_train.set_defaults(func=train)

    parser_predict = subparsers.add_parser("predict", parents=[parser], help="Module to identify antimicrobial peptides from amino acid sequence data using the random forest classifier trained using the train module.")
    requiredNamed = parser_predict.add_argument_group('required arguments')
    parser_predict.add_argument("-m", "--model", default="amPEP.model", dest="model", help="Random forest model to use for classification. Random forest model is trained and outputted using the train module. Default is amPEP.model in the current directory.")
    requiredNamed.add_argument("-i", "--input-sequences", required=True, dest="seq_file", help="Fasta file containing amino acid sequences to be classified.")
    parser_predict.add_argument("-o", "--output-file", dest="out_file", help="Ooutput file to save classification results. Default is stdout.")
    parser_predict.set_defaults(func=predict)

    if len(sys.argv) == 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        parser.print_help(sys.stderr)
        sys.exit(1)
    else:
        args = parser.parse_args()
        args.func(args)

def train(args):
    if args.max_tree <= args.min_tree:
        print(f"--max-trees ({args.max_tree}) less than or equal to --min-trees ({args.min_tree}). --max-trees must be larger than --min-trees",
              file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.positive, "r") as training_positive:
            positive_df = score(training_positive)
            positive_df['classi'] = 1
    except FileNotFoundError:
        print(f"AMP positive fasta sequence file: {args.positive} not found!", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.negative, "r") as training_negative:
            negative_df = score(training_negative)
            negative_df['classi'] = 0
    except FileNotFoundError:
        print(f"AMP negative fasta sequence file: {args.negative} not found!", file=sys.stderr)
        sys.exit(1)

    feature_drop_list = []
    if args.drop_feature:
        try:
            with open(args.drop_feature, "r") as drop_data:
                for feature in drop_data:
                    feature = feature.strip()
                    feature_drop_list.append(feature)
        except FileNotFoundError:
            print(f"Feature drop file: {args.drop_feature} not found!", file=sys.stderr)
            sys.exit(1)

    feature_drop_list.append("classi")
    training_df = pd.concat([positive_df, negative_df])
    training_df = sklearn.utils.shuffle(training_df, random_state=args.seed)
    X = training_df.drop(columns=feature_drop_list)
    #X = training_df.drop(columns=['classi', 'polarizability.1.0'])
    y = training_df.classi
    clf = RandomForestClassifier(n_estimators=args.num_trees, oob_score=True,
                                 random_state=args.seed,
                                 n_jobs=args.num_processes)
    clf.fit(X, y)
    print(f"Out-of-bag accuracy: {clf.oob_score_}", file=sys.stderr)
    #pred_train = np.argmax(clf.oob_decision_function_, axis=1).tolist()
    #train_name = training_df.index.tolist()
    #rates_df = pd.DataFrame(list(zip(pred_train, train_name)), columns=["prediction", "name"])
    #rates_df.to_csv("oob_classify_random.csv", index=False)
    #print(metrics.roc_auc_score(y, pred_train))
    if args.tree_test:
        min_estimators = args.min_tree
        max_estimators = args.max_tree
        print("n_estimators\toob_error")
        for i in range(min_estimators, max_estimators + 1):
            clf_ = RandomForestClassifier(n_estimators=i, oob_score=True,
                                          random_state=args.seed,
                                          n_jobs=args.num_processes)
            clf_.fit(X, y)
            oob_error = 1 - clf_.oob_score_
            print(f"{i}\t{oob_error}")

    if args.feature_importance:
        oob_dropcol_importances(clf, X, y, args.seed, args.num_processes)

    try:
        with open("amPEP.model", "wb") as model_pickle:
            pickle.dump(clf, model_pickle)
    except IOError:
        print("Error in writing model to file!", file=sys.stderr)

def predict(args):
    try:
        with open(args.model, "rb") as model_handle:
            clf = pickle.load(model_handle)
    except FileNotFoundError:
        print(f"Model file: {args.model} not found!", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.seq_file, "r") as classify_input:
            classify_df = score(classify_input)
    except FileNotFoundError:
        print(f"Sequence file: {args.seq_file} not found!", file=sys.stderr)
        sys.exit(1)

    if args.out_file:
        classify_output = open(args.out_file, "w")

    feature_drop_list = []
    if args.drop_feature:
        try:
            with open(args.drop_feature, "r") as drop_data:
                for feature in drop_data:
                    feature = feature.strip()
                    feature_drop_list.append(feature)
        except FileNotFoundError:
            print(f"Feature drop file: {args.drop_feature} not found!", file=sys.stderr)
            sys.exit(1)

    if feature_drop_list:
        classify_df = classify_df.drop(columns=feature_drop_list)
    id_info = classify_df.index.tolist()

    if args.out_file:
        print("probability_nonAMP\tprobability_AMP\tpredicted\tseq_id", file=classify_output)
    else:
        print("probability_nonAMP\tprobability_AMP\tpredicted\tseq_id")
    preds = clf.predict_proba(classify_df)
    for i, pred in enumerate(preds):
        pred_list = pred.tolist()
        if clf.predict(classify_df.loc[id_info[i], :].to_numpy().reshape(1, -1))[0] == 1:
            predicted = "AMP"
        else:
            predicted = "nonAMP"
        output_line = "{}\t{}\t{}".format("\t".join([str(y) for y in pred_list]),
                                          predicted,
                                          id_info[i])
        if args.out_file:
            print(output_line, file=classify_output)
        else:
            print(output_line)
    if args.out_file:
        classify_output.close()

def score(fasta_handle):
    CTD = {'hydrophobicity': {1: ['R', 'K', 'E', 'D', 'Q', 'N'], 2: ['G', 'A', 'S', 'T', 'P', 'H', 'Y'], 3: ['C', 'L', 'V', 'I', 'M', 'F', 'W']},
           'normalized.van.der.waals': {1: ['G', 'A', 'S', 'T', 'P', 'D', 'C'], 2: ['N', 'V', 'E', 'Q', 'I', 'L'], 3: ['M', 'H', 'K', 'F', 'R', 'Y', 'W']},
           'polarity': {1: ['L', 'I', 'F', 'W', 'C', 'M', 'V', 'Y'], 2: ['P', 'A', 'T', 'G', 'S'], 3: ['H', 'Q', 'R', 'K', 'N', 'E', 'D']},
           'polarizability': {1: ['G', 'A', 'S', 'D', 'T'], 2: ['C', 'P', 'N', 'V', 'E', 'Q', 'I', 'L'], 3: ['K', 'M', 'H', 'F', 'R', 'Y', 'W']},
           'charge': {1: ['K', 'R'], 2: ['A', 'N', 'C', 'Q', 'G', 'H', 'I', 'L', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'], 3: ['D', 'E']},
           'secondary': {1: ['E', 'A', 'L', 'M', 'Q', 'K', 'R', 'H'], 2: ['V', 'I', 'Y', 'C', 'W', 'F', 'T'], 3: ['G', 'N', 'P', 'S', 'D']},
           'solvent': {1: ['A', 'L', 'F', 'C', 'G', 'I', 'V', 'W'], 2: ['R', 'K', 'Q', 'E', 'N', 'D'], 3: ['M', 'S', 'P', 'T', 'H', 'Y']}}
    header = []
    groups = [1, 2, 3]
    values = [0, 25, 50, 75, 100]
    for AAproperty in CTD:
        for types in groups:
            for numbers in values:
                label = ""
                label = label.join("{}.{}.{}".format(AAproperty, types, numbers))
                header.append(label)
    All_groups = []
    Sequence_names = []
    for sequences in SeqIO.parse(fasta_handle, "fasta"):
        sequence_name = sequences.id
        Sequence_names.append(sequence_name)
        sequence = str(sequences.seq)
        sequencelength = len(sequence)
        Sequence_group = []
        for AAproperty in CTD:
            propvalues = ""
            for letter in sequence:
                if letter in CTD[AAproperty][1]:
                    propvalues += "1"
                elif letter in CTD[AAproperty][2]:
                    propvalues += "2"
                else:
                    propvalues += "3"
            abpos_1 = [i for i in range(len(propvalues)) if propvalues.startswith("1", i)]
            abpos_1 = [x+1 for x in abpos_1]
            abpos_1.insert(0, "-")
            abpos_2 = [i for i in range(len(propvalues)) if propvalues.startswith("2", i)]
            abpos_2 = [x+1 for x in abpos_2]
            abpos_2.insert(0, "-")
            abpos_3 = [i for i in range(len(propvalues)) if propvalues.startswith("3", i)]
            abpos_3 = [x+1 for x in abpos_3]
            abpos_3.insert(0, "-")
            property_group1_length = propvalues.count("1")
            if property_group1_length == 0:
                Sequence_group.extend([0, 0, 0, 0, 0])
            elif property_group1_length == 1:
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
            elif property_group1_length == 2:
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[round((0.5*property_group1_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_1[round((0.75*property_group1_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_1[property_group1_length]/sequencelength)*100)
            else:
                Sequence_group.append((abpos_1[1]/sequencelength)*100)
                Sequence_group.append((abpos_1[round((0.25*property_group1_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_1[round((0.5*property_group1_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_1[round((0.75*property_group1_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_1[property_group1_length]/sequencelength)*100)

            property_group2_length = propvalues.count("2")
            if property_group2_length == 0:
                Sequence_group.extend([0, 0, 0, 0, 0])
            elif property_group2_length == 1:
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
            elif property_group2_length == 2:
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[round((0.5*property_group2_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_2[round((0.75*property_group2_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_2[property_group2_length]/sequencelength)*100)
            else:
                Sequence_group.append((abpos_2[1]/sequencelength)*100)
                Sequence_group.append((abpos_2[round((0.25*property_group2_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_2[round((0.5*property_group2_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_2[round((0.75*property_group2_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_2[property_group2_length]/sequencelength)*100)

            property_group3_length = propvalues.count("3")
            if property_group3_length == 0:
                Sequence_group.extend([0, 0, 0, 0, 0])
            elif property_group3_length == 1:
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
            elif property_group3_length == 2:
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[round((0.5*property_group3_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_3[round((0.75*property_group3_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_3[property_group3_length]/sequencelength)*100)
            else:
                Sequence_group.append((abpos_3[1]/sequencelength)*100)
                Sequence_group.append((abpos_3[round((0.25*property_group3_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_3[round((0.5*property_group3_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_3[round((0.75*property_group3_length)-0.1)]/sequencelength)*100)
                Sequence_group.append((abpos_3[property_group3_length]/sequencelength)*100)
        All_groups.append(Sequence_group)

    Property_dataframe = pd.DataFrame.from_dict(All_groups)
    Property_dataframe.columns = header
    Property_dataframe.index = Sequence_names
    return Property_dataframe

def oob_dropcol_importances(rf, X_train, y_train, seed, n_jobs):
    """
    Compute drop-column feature importances for scikit-learn.
    Given a RandomForestClassifier or RandomForestRegressor in rf
    and training X and y data, return a data frame with columns
    Feature and Importance sorted in reverse order by importance.
    A clone of rf is trained once to get the baseline score and then
    again, once per feature to compute the drop in out of bag (OOB)
    score.
    return: A data frame with Feature, Importance columns
    SAMPLE CODE
    rf = RandomForestRegressor(n_estimators=100, n_jobs=-1, oob_score=True)
    X_train, y_train = ..., ...
    rf.fit(X_train, y_train)
    imp = oob_dropcol_importances(rf, X_train, y_train)
    """
    rf_ = clone(rf)
    rf_.random_state = seed
    rf_.oob_score = True
    rf_.n_jobs = n_jobs
    rf_.fit(X_train, y_train)
    baseline = rf_.oob_score_
    imp = []
    for col in X_train.columns:
        rf_ = clone(rf)
        rf_.random_state = seed
        rf_.oob_score = True
        rf_.n_jobs = n_jobs
        rf_.fit(X_train.drop(col, axis=1), y_train)
        drop_in_score = baseline - rf_.oob_score_
        imp.append(drop_in_score)
    imp = np.array(imp)
    I = pd.DataFrame(data={'Feature':X_train.columns, 'Importance':imp})
    I = I.set_index('Feature')
    I = I.sort_values('Importance', ascending=False)
    I.to_csv("feature.importances.dropcol.oob.csv", index=True)

if __name__ == "__main__":
    main()
