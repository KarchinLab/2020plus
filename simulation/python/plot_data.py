import utils.python.plot as myplt

def oncogene_precision_errorbar(panel_dict, save_path,
                                title='', xlabel='', ylabel=''):
    path_holder = ''
    counter, dict_len,  = 0, len(panel_dict)
    for label, wp in panel_dict.iteritems():
        samprates = wp.items
        mean_prec = wp.ix[:, 'oncogene', 'precision mean']
        sem_prec = wp.ix[:, 'oncogene', 'precision sem']
        ci_prec = 1.96 * sem_prec  # confidence interval

        # figure out whether to save yet
        counter += 1  # counter for figuring out when to save plot
        if counter == dict_len:
            path_holder = save_path

        myplt.errorbars(samprates, mean_prec, ci_prec,
                        save_path=path_holder,
                        title=title, xlabel=xlabel,
                        ylabel=ylabel, label=label)


def oncogene_recall_errorbar(panel_dict, save_path,
                             title='', xlabel='', ylabel=''):
    path_holder = ''
    counter, dict_len,  = 0, len(panel_dict)
    for label, wp in panel_dict.iteritems():
        samprates = wp.items
        mean_rec = wp.ix[:, 'oncogene', 'recall mean']
        sem_rec = wp.ix[:, 'oncogene', 'recall sem']
        ci_rec = 1.96 * sem_rec  # confidence interval

        # figure out whether to save yet
        counter += 1  # counter for figuring out when to save plot
        if counter == dict_len:
            path_holder = save_path

        myplt.errorbars(samprates, mean_rec, ci_rec,
                        save_path=path_holder,
                        title=title, xlabel=xlabel,
                        ylabel=ylabel, label=label)


def oncogene_pr_auc_errorbar(panel_dict, save_path,
                             title='', xlabel='', ylabel=''):
    path_holder = ''
    counter, dict_len,  = 0, len(panel_dict)
    for label, wp in panel_dict.iteritems():
        samprates = wp.items
        mean_rec = wp.ix[:, 'oncogene', 'PR AUC mean']
        sem_rec = wp.ix[:, 'oncogene', 'PR AUC sem']
        ci_rec = 1.96 * sem_rec  # confidence interval

        # figure out whether to save yet
        counter += 1  # counter for figuring out when to save plot
        if counter == dict_len:
            path_holder = save_path

        myplt.errorbars(samprates, mean_rec, ci_rec,
                        save_path=path_holder,
                        title=title, xlabel=xlabel,
                        ylabel=ylabel, label=label)


def oncogene_roc_auc_errorbar(panel_dict, save_path,
                              title='', xlabel='', ylabel=''):
    path_holder = ''
    counter, dict_len,  = 0, len(panel_dict)
    for label, wp in panel_dict.iteritems():
        samprates = wp.items
        mean_rec = wp.ix[:, 'oncogene', 'ROC AUC mean']
        sem_rec = wp.ix[:, 'oncogene', 'ROC AUC sem']
        ci_rec = 1.96 * sem_rec  # confidence interval

        # figure out whether to save yet
        counter += 1  # counter for figuring out when to save plot
        if counter == dict_len:
            path_holder = save_path

        myplt.errorbars(samprates, mean_rec, ci_rec,
                        save_path=path_holder,
                        title=title, xlabel=xlabel,
                        ylabel=ylabel, label=label)
