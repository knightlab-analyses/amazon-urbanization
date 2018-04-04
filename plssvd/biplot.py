from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import matplotlib.cm as cmx
import numpy as np
import pandas as pd
from collections import OrderedDict


def make_biplot(samples,
                features=None,
                sample_metadata=None,
                feature_metadata=None,
                sample_color_category=None,
                feature_color_category=None,
                sample_color_dict=None,
                feature_color_dict=None,
                sample_zorder=None,
                feature_zorder=None,
                sample_labels=None,
                feature_labels=None,
                **kwargs):

    figure_size = (25, 25)
    samples_x = 'PCA1'
    samples_y = 'PCA2'
    samp_col = 'RdGy'
    samp_alpha = 1
    samp_marker = 'o'
    samp_ms = 8
    samp_leg_loc = 2
    features_x = 'PCA1'
    features_y = 'PCA2'
    feat_col = 'Set1'
    feat_alpha = 1
    arrow_width = 0.02
    arrow_head = 0.05
    feat_leg_loc = 1
    feature_order = 0
    sample_drop_list = []
    show_color_drop = False
    sample_drop_col = ['#FFFFFF']
    eigenvalues = []
    x_pad = 0.3
    y_pad = 0.3

    for key, value in kwargs.items():
        if key == 'figure_size':
            figure_size = value
        if key == 'samples_x':
            samples_x = value
        if key == 'samples_y':
            samples_y = value
        if key == 'samp_col':
            samp_col = value
        if key == 'samp_alpha':
            samp_alpha = value
        if key == 'samp_marker':
            samp_marker = value
        if key == 'samp_ms':
            samp_ms = value
        if key == 'samp_leg_loc':
            samp_leg_loc = value
        if key == 'features_x':
            samples_x = value
        if key == 'features_y':
            samples_y = value
        if key == 'feat_col':
            feat_col = value
        if key == 'feat_alpha':
            feat_alpha = value
        if key == 'arrow_width':
            arrow_width = value
        if key == 'arrow_head':
            arrow_head = value
        if key == 'feat_leg_loc':
            feat_leg_loc = value
        if key == 'feature_order':
            if value == 0:
                feature_order = 0
            if value == 1:
                feature_order = 1
        if key == 'sample_drop_list':
            sample_drop_list = value
        if key == 'show_color_drop':
            show_color_drop = value
        if key == 'sample_drop_col':
            sample_drop_col = value
        if key == 'eigenvalues':
            eigenvalues = value
        if key == 'x_pad':
            x_pad = value
        if key == 'y_pad':
            y_pad = value

    if not isinstance(samples, pd.core.frame.DataFrame):
        raise ValueError('`samples` must be a `pd.DataFrame`, '
                         'not %r.' % type(samples).__name__)

    if features is not None:
        if not isinstance(features, pd.core.frame.DataFrame):
            raise ValueError('`features` must be a `pd.DataFrame`, '
                             'not %r.' % type(features).__name__)

    if sample_metadata is not None:
        if not isinstance(sample_metadata, pd.core.frame.DataFrame):
            raise ValueError('`sample_metadata` must be a `pd.DataFrame`, '
                             'not %r.' % type(sample_metadata).__name__)

    if feature_metadata is not None:
        if not isinstance(feature_metadata, pd.core.frame.DataFrame):
            raise ValueError('`feature_metadata` must be a `pd.DataFrame`, '
                             'not %r.' % type(feature_metadata).__name__)

    if sample_color_dict is not None:
        if not isinstance(sample_color_dict, dict):
            raise ValueError('`sample_color_dict` must be a `dictionary`, '
                             'not %r.' % type(sample_color_dict).__name__)

    if feature_color_dict is not None:
        if not isinstance(feature_color_dict, dict):
            raise ValueError('`feature_color_dict` must be a `dictionary`, '
                             'not %r.' % type(feature_color_dict).__name__)

    if sample_metadata is not None and sample_color_dict is None:
        if sample_color_category is None:
            raise ValueError('sample_color_category must be a specified')

    if sample_metadata is not None and sample_color_dict is not None:
        if sample_color_category is None:
            raise ValueError('sample_color_category must be a specified')

    if feature_metadata is not None and feature_color_dict is not None:
        if feature_color_category is None:
            raise ValueError('feature_color_category must be a specified')

    if sample_drop_list is not None:
        if not isinstance(sample_drop_list, list):
            raise ValueError('`sample_drop_list` must be a `list`, '
                             'not %r.' % type(sample_drop_list).__name__)

    if sample_drop_col is not None:
        if not isinstance(sample_drop_col, list):
            raise ValueError('`sample_drop_col` must be a `list`, '
                             'not %r.' % type(sample_drop_col).__name__)

    if sample_metadata is not None:
        if (samples.index != sample_metadata.index).any():
            samples = samples.sort_index(axis=0)
            sample_metadata = sample_metadata.sort_index(axis=0)

    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111)

    sample_colors = plt.get_cmap(samp_col)
    feature_colors = plt.get_cmap(feat_col)

    sample_group_append = []
    colorVal = []

    if sample_metadata is None:
        ax.plot(np.ravel(samples[samples_x]),
                np.ravel(samples[samples_y]),
                marker=samp_marker, linestyle='',
                ms=samp_ms, alpha=samp_alpha)

    if (sample_metadata is not None and sample_color_dict is None):
        sample_groups = samples.groupby(sample_metadata[sample_color_category])

        if len(sample_drop_list) > 0:

            def dropf(x):
                return x not in sample_drop_list
            index_drop = sample_metadata[sample_color_category].apply(dropf)

            samp_r = samples.loc[sample_metadata.index]
            samp_met_r = sample_metadata.loc[index_drop][sample_color_category]

            for name, group in samp_r.groupby(samp_met_r):
                sample_group_append.append(name)
                sample_group_append = sorted(list(set(sample_group_append)))
                cNorm = colors.Normalize(vmin=0,
                                         vmax=(len(sample_group_append)-1))
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=sample_colors)

            for index, row in enumerate(sample_group_append):
                colorVal.append(scalarMap.to_rgba(index))

            if not show_color_drop:
                sample_color_dict = dict(zip(sample_group_append, colorVal))
                sample_color_dict = OrderedDict(
                    sorted(sample_color_dict.items(),
                           key=lambda x: x[0],
                           reverse=True))

                for name, group in samp_r.groupby(samp_met_r):
                    kwds = {}
                    if sample_zorder is not None:
                        kwds['zorder'] = sample_zorder[name]
                    ax.plot(np.ravel(group[samples_x]),
                            np.ravel(group[samples_y]),
                            marker=samp_marker, linestyle='', ms=samp_ms,
                            color=sample_color_dict[name],
                            label=name, alpha=samp_alpha,
                            **kwds)

            else:
                color_drop_append = []
                if len(sample_drop_col) == 1:
                    for index in range(len(sample_drop_list)):
                        color_drop_append.append(sample_drop_col[0])
                    colorVal = colorVal + color_drop_append

                if len(sample_drop_col) == len(sample_drop_list):
                    for index in range(len(sample_drop_list)):
                        color_drop_append.append(sample_drop_col[index])
                    colorVal = colorVal + color_drop_append

                sample_group_append = list(sample_group_append)
                sample_group_append += list(sample_drop_list)
                sample_color_dict = dict(zip(sample_group_append, colorVal))
                sample_color_dict = OrderedDict(
                    sorted(sample_color_dict.items(),
                           key=lambda x: x[0],
                           reverse=True))

                for name, group in sample_groups:
                    kwds = {}
                    if sample_zorder is not None:
                        kwds['zorder'] = sample_zorder[name]
                    if name not in sample_drop_list:
                        ax.plot(np.ravel(group[samples_x]),
                                np.ravel(group[samples_y]),
                                marker=samp_marker, linestyle='', ms=samp_ms,
                                color=sample_color_dict[name],
                                label=name, alpha=samp_alpha,
                                **kwds)

                for name, group in sample_groups:
                    kwds = {}
                    if sample_zorder is not None:
                        kwds['zorder'] = sample_zorder[name]
                    if name in sample_drop_list:
                        ax.plot(np.ravel(group[samples_x]),
                                np.ravel(group[samples_y]),
                                marker='samp_marker', linestyle='', ms=samp_ms,
                                color=sample_color_dict[name],
                                label=name, alpha=samp_alpha,
                                **kwds)

        else:
            sample_group_append = []
            for name, group in sample_groups:
                sample_group_append.append(name)
                sample_group_append = sorted(list(set(sample_group_append)))
                cNorm = colors.Normalize(vmin=0,
                                         vmax=(len(sample_group_append)-1))
                scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=sample_colors)

            for index, row in enumerate(sample_group_append):
                colorVal.append(scalarMap.to_rgba(index))

            sample_color_dict = dict(zip(sample_group_append, colorVal))
            sample_color_dict = OrderedDict(
                sorted(sample_color_dict.items(),
                       key=lambda x: x[0],
                       reverse=True))

            for name, group in sample_groups:
                kwds = {}
                if sample_zorder is not None:
                    kwds['zorder'] = sample_zorder[name]
                ax.plot(np.ravel(group[samples_x]),
                        np.ravel(group[samples_y]),
                        marker=samp_marker, linestyle='', ms=samp_ms,
                        color=sample_color_dict[name],
                        label=name, alpha=samp_alpha,
                        **kwds)

        sample_color_dict = None

    if (sample_metadata is not None and sample_color_dict is not None):
        sample_groups = samples.groupby(sample_metadata[sample_color_category])

        if len(sample_drop_list) > 0:

            def dropf(x):
                return x not in sample_drop_list
            index_drop = sample_metadata[sample_color_category].apply(dropf)

            samp_r = samples.loc[sample_metadata.index]
            samp_met_r = sample_metadata.loc[index_drop]

            sample_color_dict = OrderedDict(
                sorted(sample_color_dict.items(),
                       key=lambda x: x[0],
                       reverse=True))

            sample_groups = samp_r.groupby(samp_met_r[sample_color_category])
            for name, group in sample_groups:
                kwds = {}
                if sample_zorder is not None:
                    kwds['zorder'] = sample_zorder[name]
                ax.plot(np.ravel(group[samples_x]),
                        np.ravel(group[samples_y]),
                        marker=samp_marker, linestyle='', ms=samp_ms,
                        color=sample_color_dict[name],
                        label=name, alpha=samp_alpha,
                        **kwds)

        if not sample_drop_list:
            sample_color_dict = OrderedDict(
                sorted(sample_color_dict.items(),
                       key=lambda x: x[0],
                       reverse=True))

            for name, group in sample_groups:
                kwds = {}
                if sample_zorder is not None:
                    kwds['zorder'] = sample_zorder[name]
                ax.plot(np.ravel(group[samples_x]),
                        np.ravel(group[samples_y]),
                        marker=samp_marker, linestyle='', ms=samp_ms,
                        color=sample_color_dict[name],
                        label=name, alpha=samp_alpha,
                        **kwds)
                for i in range(group[features_x].shape[0]):
                    if name != 'None':
                        _id = group.index[i]
                        _label = (sample_labels[_id] if _id in sample_labels
                                  else 'Unknown') if sample_labels else _id
                        ax.text(group.loc[_id, features_x],
                                group.loc[_id, features_y], _label)

        sample_color_dict = None

    ax2 = ax.twinx()

    if sample_color_category is not None:
        ax.legend(title=sample_color_category, loc=samp_leg_loc, numpoints=1)
    else:
        ax.legend(loc=samp_leg_loc, numpoints=1)

    ax2.set_ylim(ax.get_ylim())

    recs = []
    feature = []
    otu_feature_append = []
    colorVal = []

    if (features is not None and feature_metadata is None):
        for index, row in features.iterrows():
            ax2.arrow(0, 0, row[features_x], row[features_y],
                      width=arrow_width, head_width=arrow_head,
                      alpha=feat_alpha, color='b')

    if (features is not None and
            feature_metadata is not None and
            feature_color_category is None):
            otu_feature_append = []

            feature_groups = features.groupby(feature_metadata.columns[0])
            for name, group in feature_groups:
                otu_feature_append.append(name)

            otu_feature_append = sorted(list(set(otu_feature_append)))
            cNorm = colors.Normalize(vmin=0, vmax=(len(otu_feature_append)-1))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=feature_colors)

            for index, row in enumerate(otu_feature_append):
                colorVal.append(scalarMap.to_rgba(index))

            feature_color_dict = dict(zip(otu_feature_append, colorVal))
            feature_color_dict = OrderedDict(
                sorted(feature_color_dict.items(),
                       key=lambda x: x[0]))

            for name, group in feature_groups:
                for i in range(group[features_x].shape[0]):
                    _id = group.index[i]
                    ax2.arrow(0, 0,
                              group.loc[_id, features_x],
                              group.loc[_id, features_y],
                              width=arrow_width, head_width=arrow_head,
                              alpha=feat_alpha,
                              color=feature_color_dict[name])
                    if name != 'None':
                        _label = (feature_labels[_id] if _id in feature_labels
                                  else 'Unknown') if feature_labels else _id
                        ax2.text(group.loc[_id, features_x],
                                 group.loc[_id, features_y],
                                 _label)

            for key in (sorted(feature_zorder, key=feature_zorder.get)
                        if feature_zorder else feature_color_dict):
                recs.append(mpatches.Rectangle((0, 0), 1, 1,
                            fc=feature_color_dict[key],
                            alpha=feat_alpha))
                feature.append(key)

            ax2.legend(recs, feature, loc=feat_leg_loc,
                       title=feature_color_category)
            feature_color_dict = None

    if (features is not None and
            feature_metadata is not None and
            feature_color_category is not None):

        feature_groups = features.groupby(
            feature_metadata[feature_color_category])

        if feature_color_dict is None:
            otu_feature_append = []
            feature_groups = features.groupby(
                feature_metadata[feature_color_category])

            for name, group in feature_groups:
                otu_feature_append.append(name)

            otu_feature_append = sorted(list(set(otu_feature_append)))
            cNorm = colors.Normalize(vmin=0,
                                     vmax=(len(otu_feature_append)-1))
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=feature_colors)

            for index, row in enumerate(otu_feature_append):
                colorVal.append(scalarMap.to_rgba(index))
            feature_color_dict = dict(zip(otu_feature_append, colorVal))
            feature_color_dict = OrderedDict(
                sorted(feature_color_dict.items(),
                       key=lambda x: x[0]))

        for name, group in feature_groups:
            for i in range(group[features_x].shape[0]):
                _id = group.index[i]

                kwds = {}
                if feature_zorder is not None:
                    kwds['zorder'] = feature_zorder[name]
                ax2.arrow(0, 0,
                          group.loc[_id, features_x],
                          group.loc[_id, features_y],
                          width=arrow_width, head_width=arrow_head,
                          alpha=feat_alpha,
                          color=feature_color_dict[name],
                          **kwds)
                if name != 'None':
                    _label = (feature_labels[_id] if _id in feature_labels
                              else 'Unknown') if feature_labels else _id
                    ax2.text(group.loc[_id, features_x],
                             group.loc[_id, features_y],
                             _label)

        for key in (sorted(feature_zorder, key=feature_zorder.get)
                    if feature_zorder else feature_color_dict):
            recs.append(mpatches.Rectangle((0, 0), 1, 1,
                                           fc=feature_color_dict[key],
                                           alpha=feat_alpha))
            feature.append(key)
        ax2.legend(recs, feature, loc=feat_leg_loc,
                   title=feature_color_category)

    if features is not None:
        xmin = min([min(samples.ix[:, 0]), min(features.ix[:, 0])])
        xmax = max([max(samples.ix[:, 0]), max(features.ix[:, 0])])
        ymin = min([min(samples.ix[:, 1]), min(features.ix[:, 1])])
        ymax = max([max(samples.ix[:, 1]), max(features.ix[:, 1])])
        xpad = (xmax - xmin) * x_pad
        ypad = (ymax - ymin) * y_pad

        ax.set_zorder(ax2.get_zorder()+(1-feature_order))
        ax.patch.set_visible(False)

        ax.set_xlim(xmin - xpad, xmax + xpad)
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ax2.set_xlim(xmin - xpad, xmax + xpad)
        ax2.set_ylim(ymin - ypad, ymax + ypad)
        ax2.set_yticks([])
    else:
        xmin = min([min(samples.ix[:, 0])])
        xmax = max([max(samples.ix[:, 0])])
        ymin = min([min(samples.ix[:, 1])])
        ymax = max([max(samples.ix[:, 1])])
        xpad = (xmax - xmin) * x_pad
        ypad = (ymax - ymin) * y_pad

        ax.set_xlim(xmin - xpad, xmax + xpad)
        ax.set_ylim(ymin - ypad, ymax + ypad)
        ax2.set_yticks([])
    if len(eigenvalues) > 2:
        e_0 = eigenvalues[0]
        e_1 = eigenvalues[1]
        ax.set_xlabel('PC 1 ({:.2%})'.format(e_0**2/sum(eigenvalues**2)))
        ax.set_ylabel('PC 2 ({:.2%})'.format(e_1**2/sum(eigenvalues**2)))

    return fig, [ax, ax2]
