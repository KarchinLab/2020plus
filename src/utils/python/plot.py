"""
The plot.py does the actual calling of plot commands from matplotlib.
Essentially, plot.py encapsulates all the minor tweaks needed in matplotlib
to make a reasonable looking plot.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import IPython


def heatmap(df, file_path, title='', xlabel='', ylabel='', cmap=plt.cm.Blues):
    """Plot a heatmap from a pandas dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        data for heatmap plotting
    file_path : str
        path to save figure (png, pdf, etc.)
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    cmap : cm
        color scheme for heatmap

    Code from:
    http://stackoverflow.com/questions/14391959/heatmap-in-matplotlib-with-pcolor

    For information on colorbars:
    http://stackoverflow.com/questions/13943217/how-to-add-colorbars-to-scatterplots-created-like-this
    """
    df = df.fillna(0)  # fills missing values with 0's

    # make heatmap
    fig, ax = plt.subplots()
    hmap = ax.pcolor(df, cmap=cmap, alpha=0.8)

    fig = plt.gcf()
    fig.set_size_inches(8,11)

    # turn off the frame
    ax.set_frame_on(False)

    # put the major ticks at the middle of each cell
    ax.set_yticks(np.arange(df.shape[0])+0.5, minor=False)
    ax.set_xticks(np.arange(df.shape[1])+0.5, minor=False)

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    # Set the labels
    labels = df.index
    ax.set_xticklabels(labels, minor=False)
    ax.set_yticklabels(df.index, minor=False)

    # rotate the
    plt.xticks(rotation=90)

    ax.grid(False)

    # Turn off all the ticks
    ax = plt.gca()

    for t in ax.xaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False
    for t in ax.yaxis.get_major_ticks():
        t.tick1On = False
        t.tick2On = False

    # handle labels
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()

    # turn on color bar
    cb = plt.colorbar(hmap)
    cb.set_label('Transition Probability')

    # save figure
    plt.savefig(file_path)
    plt.close()


def barplot(df,
            file_path,
            kind='bar',
            yerr=None,
            xerr=None,
            ecolor=None,
            title='',
            xlabel='',
            ylabel='',
            stacked=False):
    """barplot generates/saves a bar plot from a pandas data frame.

    Parameters
    ----------
    df : pd.DataFrame
        data frame for bar plot
    file_path : str
        path to save bar plot figure
    kind : str, ['bar' | 'barh']
        vertical or horizontal bar chart
    stderror : list
        stdev of each bar
    Matplotlib options for plotting
    """
    if yerr is not None:
        # error bars for y-axis
        df.plot(kind=kind, yerr=yerr, ecolor=ecolor, stacked=stacked)
    elif xerr is not None:
        # error bars for x-axis
        df.plot(kind=kind, xerr=xerr, ecolor=ecolor, stacked=stacked)
    else:
        # normal bar plot
        df.plot(kind=kind, stacked=stacked)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(file_path)
    plt.clf()  # clear figure
    plt.close()


def histogram(df,
              save_path,
              bins=False,
              log=False,
              title='',
              xlabel='',
              ylabel=''):
    """Plots a histogram using matplotlib.

    Parameters
    ----------
    df : pd.DataFrame
        one dimensional data frame or series
    file_path : str
        path to save figure
    bins : list
        bin positions for histogram
    log : Bool
        boolean for log scaling y-axis
    title : str
        title of plot
    xlabel : str
        label on x-axis
    ylabel : str
        label on y-axis
    """
    if bins:
        df.hist(bins=bins, log=log)
    else:
        df.hist(log=log)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.clf()  # clear figure
    plt.close()


def line(data,
         file_path,
         style=[],
         title='',
         xlabel='',
         ylabel='',
         logx=False,
         logy=False,
         vlines=[]):
    """Plots a line plot using matplotlib.

    Parameters
    ----------
    data : pd.DataFrame
        two column df with x and y values
    file_path : str
        path to save figure
    title : str
        graph title
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    vlines : list of int
        draw vertical lines at positions
    """
    # plot data
    data.plot(kind='line', style=style)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    # log scale if neccessary
    if logx:
        plt.xscale('log')
    if logy:
        plt.yscale('log')

    # plot vertical lines
    ymin, ymax = plt.ylim()  # get plotting range of y-axis
    for l in vlines:
        plt.vlines(l, ymin=ymin,
                   ymax=ymax,
                   color='red')

    plt.tight_layout()  # adjust plot margins
    plt.savefig(file_path)  # save figure
    plt.clf()  # clear figure
    plt.close()


def scatter(x, y,
            file_path,
            colors='',
            size=20,
            title='',
            xlabel='',
            ylabel=''):
    """Create a 2D scatter plot. Many of the optional arguements
    deal with formatting the plot.

    Parameters
    ----------
    x : list|array
        container for x-axis data
    y : list|array
        container for y-axis data
    file_path : str
        path to save figure
    colors : str|list, default=''
        either single color (e.g. 'blue') or a list
    size : int|list, default=20
        int for marker size or a list of ints
    title : str, default=''
        title for plot
    xlabel : str, dfault=''
        x-axis label
    ylabel : str, default=''
        y-axis label
    """
    if colors:
        # plot user's color if supplied
        plt.scatter(x, y, c=colors, s=size)
    else:
        # otherwise rely on default color
        plt.scatter(x, y, s=size)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(file_path)
    plt.clf()  # clear figure
    plt.close()


def line_fill_between(data, sem,
                      save_path,
                      style=[],
                      title='',
                      xlabel='',
                      ylabel=''):
    """Plot a line graph from data but also add a fill between effect
    base on the standard error of the mean (sem).

    Parameters
    ----------
    data : pd.DataFrame
        data for one to multiple lines
    sem : pd.DataFrame
        sem for each point in data
    save_path : str
        file path for saving
    title :str
        title of figure
    xlabel : str
        x-axis label
    ylabel : str
        y-axis label
    """
    # plot lines
    ax = data.plot(kind='line', style=style)

    # plot fill between which indicates the standard error of the mean
    line_colors = [x.get_color() for x in ax.get_lines()]
    x = data.index.values.astype(float)  # x values for plotting
    for i in range(len(line_colors)):
        y = data.iloc[:, i]
        color = line_colors[i]  # get matching line color
        tmp_sem = sem.iloc[:, i]  # get a single column
        plt.fill_between(x, y-tmp_sem, y+tmp_sem, alpha=.5, facecolor=color)

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.clf()
    plt.close()


def errorbars(x, y, err,
              save_path='',
              title='',
              xlabel='',
              ylabel='',
              label=''):
    if label:
        plt.errorbar(x, y, yerr=err, label=label, fmt='-o')
        plt.legend(loc='best')
    else:
        plt.errorbar(x, y, yerr=err, fmt='-o')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    if save_path:
        plt.savefig(save_path)
        plt.close()


def correlation_plot(x, y,
                     save_path,
                     title,
                     xlabel, ylabel):
    plt.scatter(x, y)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    line_x = np.arange(x.min(), x.max())
    line_y = slope*line_x + intercept
    plt.plot(line_x, line_y,
             label='$%.2fx + %.2f$, $R^2=%.2f$' % (slope, intercept, r_value**2))
    plt.legend(loc='best')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.clf()  # clear figure
    plt.close()


def boxplot(df, by,
            column,
            save_path,
            xlabel,
            ylabel,
            title):
    # make box plot
    if column:
        axes = df.boxplot(column=column, by=by)
    else:
        axes = df.boxplot(by=by)

    if type(axes) is np.ndarray:
        # multiple boxplot case
        multi_box_flag = True
    else:
        # only one facet for boxplot
        multi_box_flag = False

    # format plot
    if multi_box_flag:
        # plot with multiple box plot facets
        for ax in axes:
            ax.set_xlabel(xlabel)
        axes[0].set_ylabel(ylabel)
        fig = axes[0].get_figure()
        fig.suptitle('')
        plt.tight_layout()
        plt.savefig(save_path)
        plt.clf()
        plt.close()
    else:
        # plot when just 1 box plot facet
        fig = axes.get_figure()
        fig.suptitle('')  # remove auto title from pandas
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.tight_layout()
        plt.savefig(save_path)
        plt.clf()  # clear figure
        plt.close()
