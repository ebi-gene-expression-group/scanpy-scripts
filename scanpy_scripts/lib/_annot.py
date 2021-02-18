"""
Transfer annotation functions
"""

import os
import sys
import numpy as np
import numpy_groupies as npg
import scipy.sparse as sp
import pandas as pd
import joblib
from copy import deepcopy
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import adjusted_rand_score
from sklearn.utils.testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning


@ignore_warnings(category=ConvergenceWarning)
def LR_train(
        adata, groupby, use_rep='raw', use_hvg=False, use_pseudobulk=False, max_pass=20,
        save=None, model=None, **kwargs):
    groupby_var = adata.obs[groupby].cat.remove_unused_categories()
    Y = groupby_var.astype(str)
    if use_rep == 'raw':
        X = adata.raw.X
        features = adata.raw.var_names.values
    elif use_rep == 'X':
        X = adata.X
        features = adata.var_names.values
        if use_hvg and 'highly_variable' in adata.var.keys():
            k_hvg = adata.var['highly_variable'].values
            X = X[:, k_hvg]
            features = features[k_hvg]
    elif use_rep in adata.obsm.keys():
        X = adata.obsm[use_rep]
        features = np.array([f'V{i+1}' for i in range(X.shape[1])])
    else:
        raise KeyError(f'{use_rep}: invalid <use_rep>')

    if use_pseudobulk:
        if sp.issparse(X):
            summarised = np.zeros((Y.unique().size, X.shape[1]))
            for i, grp in enumerate(groupby_var.cat.categories):
                k_grp = np.where(Y == grp)[0]
                summarised[i] = np.mean(X[k_grp, :], axis=0)
            X = summarised
        else:
            X = npg.aggregate(groupby_var.cat.codes, X, axis=0, func='mean')
        Y = groupby_var.cat.categories.values

    lr = model if model is not None else LogisticRegression(penalty='l2', C=0.1, solver='saga', warm_start=True, n_jobs=-1, **kwargs)
    n_pass = 0
    while n_pass < max_pass:
        lr.fit(X, Y)
        n_pass += 1
        if lr.n_iter_ < 100:
            break
    if lr.n_iter_ >= 100:
        print('training of LR model failed to converge', file=sys.stderr)
    lr.features = features
    if save:
        outdir = os.path.dirname(save)
        if outdir:
            os.makedirs(outdir, exist_ok=True)
        joblib.dump(lr, save)
    return lr


def LR_predict(
        adata, model, use_rep='raw', use_pseudobulk=False, groupby=None, feature=None,
        key_added=None, truth=None, return_predict=False):
    if use_rep == 'raw':
        X = adata.raw.X
        features = adata.raw.var_names if feature is None else adata.raw.var[feature]
    elif use_rep == 'X':
        X = adata.X
        features = adata.var_names if feature is None else adata.var[feature]
    elif use_rep in adata.obsm.keys():
        X = adata.obsm[use_rep]
        features = np.array([f'V{i+1}' for i in range(X.shape[1])])
    else:
        raise KeyError(f'{use_rep}: invalid <use_rep>')
    features = pd.Series(features)

    if isinstance(model, str) and os.path.exists(model):
        lr = joblib.load(model)
    elif isinstance(model, LogisticRegression):
        lr = deepcopy(model)
    else:
        raise ValueError(f'{model}: invalid LR model')
    if getattr(lr, 'features', None) is None:
        if lr.n_features_in_ == features.size:
            lr.features = features.values
        else:
            raise ValueError(f'{model}: LR model has no feature names and unmatched size')

    k_x = features.isin(list(lr.features))
    print(f'{k_x.sum()} features used for prediction', file=sys.stderr)
    k_x_idx = np.where(k_x)[0]
    X = X[:, k_x_idx]
    features = features[k_x]

    ad_ft = pd.DataFrame(features.values, columns=['ad_features']).reset_index().rename(columns={'index': 'ad_idx'})
    lr_ft = pd.DataFrame(lr.features, columns=['lr_features']).reset_index().rename(columns={'index': 'lr_idx'})
    lr_idx = lr_ft.merge(ad_ft, left_on='lr_features', right_on='ad_features').sort_values(by='ad_idx').lr_idx.values

    lr.n_features_in_ = lr_idx.size
    lr.features = lr.features[lr_idx]
    lr.coef_ = lr.coef_[:, lr_idx]

    if use_pseudobulk:
        if not groupby or groupby not in adata.obs.columns:
            raise ValueError('missing or invalid `groupby`')
        groupby_var = adata.obs[groupby].cat.remove_unused_categories()
        summarised = np.zeros((groupby_var.cat.categories.size, X.shape[1]))
        for i, grp in enumerate(groupby_var.cat.categories):
            k_grp = np.where(groupby_var == grp)[0]
            if sp.issparse(X):
                summarised[i] = np.mean(X[k_grp, :], axis=0)
            else:
                summarised[i] = np.mean(X[k_grp, :], axis=0, keepdims=True)
        X = summarised

    ret = {}
    Y_predict = lr.predict(X)
    Y_prob = lr.predict_proba(X)
    max_Y_prob = Y_prob.max(axis=1)
    if use_pseudobulk:
        tmp_groupby = adata.obs[groupby].astype(str)
        tmp_prob = np.zeros(tmp_groupby.size)
        tmp_predict = tmp_prob.astype(str)
        for i, ct in enumerate(adata.obs[groupby].cat.categories):
            tmp_prob[tmp_groupby==ct] = max_Y_prob[i]
            tmp_predict[tmp_groupby==ct] = Y_predict[i]
        max_Y_prob = tmp_prob
        Y_predict = tmp_predict
    if key_added:
        adata.obs[key_added] = Y_predict
        adata.obs[key_added] = adata.obs[key_added].astype('category')
        adata.obs[f'{key_added}_prob'] = max_Y_prob
    ret['label'] = Y_predict
    ret['prob'] = pd.DataFrame(
        Y_prob,
        index=adata.obs[groupby].cat.categories if use_pseudobulk else adata.obs_names,
        columns=lr.classes_)

    if truth:
        Y_truth = adata.obs[truth].astype(str)
        ret['accuracy'] = (Y_predict == Y_truth).sum() / Y_predict.size
        ret['adjusted_rand_score'] = adjusted_rand_score(Y_truth, Y_predict)

    if return_predict:
        return ret


def annotate(adata, groupby, label, normalise_label=True, threshold=0.85, max_entry=None):
    """Annotate a clustering based on a matrix of values

    + groupby: the clustering to be annotated
    + label : cell-level annotation
    """
    annot_mat = pd.crosstab(adata.obs[groupby], adata.obs[label])
    group_size = annot_mat.sum(axis=1).values
    group_prop = annot_mat / group_size[:, np.newaxis]
    if normalise_label:
        label_size = group_prop.sum(axis=0).values
        label_prop = group_prop / label_size[np.newaxis, :]
    else:
        label_prop = group_prop
    v_max = label_prop.values.max(axis=1)

    annot_dict = {}
    for i, row in label_prop.iterrows():
        idx = np.where(row.values > row.values.max()*threshold)[0]
        od = np.argsort(row.values[idx])
        if max_entry is not None:
            max_entry = min(max_entry, len(idx))
            entries = row[idx[od]][0:max_entry]
        else:
            entries = row[idx[od]]
        gl = ';'.join(entries.index.values)
        if gl not in annot_dict:
            annot_dict[gl] = []
        annot_dict[gl].append(i)
    return annot_dict


# def LR_annotate(
#         adata,
#         train_label,
#         train_x,
#         use_rep='X',
#         highly_variable=True,
#         subset=None,
#         penalty='l1',
#         sparcity=0.2,
#         return_label=False,
# ):
#     """Annotate cells using logistic regression
#     """
#     if subset is not None:
#         train_label = train_label[subset]
#         train_x = train_x[subset, :]
#     from sklearn.linear_model import LogisticRegression
#     lr = LogisticRegression(penalty=penalty, C=sparcity)
#     lr.fit(train_x, train_label)
#     if use_rep == 'X':
#         if highly_variable:
#             predict_x = adata.X[:, adata.var.highly_variable]
#         else:
#             predict_x = adata.X
#     elif use_rep == 'raw':
#         if highly_variable:
#             predict_x = adata.raw.X[:, adata.var.highly_variable]
#         else:
#             predict_x = adata.raw.X
#     elif use_rep in adata.layers.keys():
#         predict_x = adata.layers[use_rep]
#     elif use_rep in adata.obsm.keys():
#         predict_x = adata.obsm[use_rep]
#     if sp.issparse(predict_x):
#         predict_x = predict_x.toarray()
#     if return_label:
#         return lr.predict(predict_x)
#     else:
#         return lr.predict_proba(predict_x)


# def annotate(adata, groupby, annot, annotation_matrix):
#     """Annotate a clustering based on a matrix of values
# 
#     + groupby: the clustering to be annotated
#     + annot : the annotationo to be aligned/transferred
#     + annotation_matrix: a matrix with shape (n,m) where n equals the number
#                          of clusters and m equals the number of annotation
#                          categories.
#     """
#     if groupby not in adata.obs.columns:
#         raise KeyError(f'"{groupby}" not found.')
#     if annot not in adata.obs.columns:
#         raise KeyError(f'"{annot}" not found.')
#     groupby_annot = f'{groupby}_annot'
#     clustering = adata.obs[groupby]
#     annotation = adata.obs[annot]
#     k_nan = annotation == 'nan'
#     cell_types = np.unique(annotation[~k_nan].values)
#     n_celltype = cell_types.shape[0]
#     n_cluster = np.unique(clustering).shape[0]
#     n_cell = clustering.shape[0]
# 
#     def dedup(names):
#         names_copy = names.copy()
#         count = {}
#         for i, name in enumerate(names):
#             n = count.get(name, 0) + 1
#             count[name] = n
#             if n > 1 or name in names[(i+1):]:
#                 names_copy[i] = f'{name}{n}'
#             else:
#                 names_copy[i] = name
#         return names_copy
# 
#     if (annotation_matrix.shape[0] == n_cluster and
#             annotation_matrix.shape[1] == n_celltype):
#         k_max = np.argmax(annotation_matrix, axis=1)
#         cluster_annot = cell_types[k_max]
#         adata.obs[groupby_annot] = adata.obs[groupby].cat.rename_categories(
#                 dedup(cluster_annot))
#     else:
#         raise ValueError(f'[annotation_matrix] not in compatible shape')
