#!/usr/bin/env python3
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc, confusion_matrix, classification_report
from sklearn.feature_selection import SelectKBest, f_classif


def build_ml_biomarker_module(deg_csv, counts_csv, metadata_csv, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    meta = pd.read_csv(metadata_csv, index_col=0)
    counts = pd.read_csv(counts_csv, sep='\t', skiprows=1, index_col=0)

    sample_cols = [c for c in counts.columns if c not in ['Chr', 'Start', 'End', 'Strand', 'Length']]
    expr_matrix = counts[sample_cols].copy()
    expr_matrix.columns = expr_matrix.columns.str.replace(r'.sortedByCoord.out.bam$', '', regex=True)

    cpm = expr_matrix.div(expr_matrix.sum(axis=0), axis=1) * 1e6
    cpm = np.log2(cpm + 1).T

    cpm = cpm.loc[meta.index]
    y = meta['condition'].map({'Normal': 0, 'Tumor': 1}).values

    selector = SelectKBest(score_func=f_classif, k=min(100, cpm.shape[1]))
    X_selected = selector.fit_transform(cpm, y)
    selected_features = cpm.columns[selector.get_support()]

    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
    classifier = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)
    oof_preds = np.zeros(len(y))

    fig, ax = plt.subplots(figsize=(6, 6))

    for i, (train, test) in enumerate(cv.split(X_selected, y)):
        classifier.fit(X_selected[train], y[train])
        viz = roc_curve(y[test], classifier.predict_proba(X_selected[test])[:, 1])
        oof_preds[test] = classifier.predict(X_selected[test])

        interp_tpr = np.interp(mean_fpr, viz[0], viz[1])
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(auc(viz[0], viz[1]))
        ax.plot(viz[0], viz[1], alpha=0.3, label=f'Fold {i+1} (AUC = {auc(viz[0], viz[1]):.2f})')

    ax.plot([0, 1], [0, 1], linestyle='--', color='r', label='Chance', alpha=0.8)
    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)

    ax.plot(mean_fpr, mean_tpr, color='b', label=f'Mean ROC (AUC = {mean_auc:.2f})', lw=2, alpha=0.8)
    ax.set(xlabel='False Positive Rate', ylabel='True Positive Rate', title='Biomarker Classifier Performance')
    ax.legend(loc="lower right")
    plt.savefig(os.path.join(out_dir, 'ml_biomarker_roc.png'), dpi=300)
    plt.close()

    classifier.fit(X_selected, y)
    importances = pd.DataFrame({
        'Gene': selected_features,
        'Importance': classifier.feature_importances_
    }).sort_values(by='Importance', ascending=False)

    importances.to_csv(os.path.join(out_dir, 'biomarker_importances.csv'), index=False)
    print("Machine learning stratification module completed successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python ml_biomarker_classifier.py <deg_csv> <counts_matrix> <metadata_csv> <out_dir>")
        sys.exit(1)
    build_ml_biomarker_module(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
