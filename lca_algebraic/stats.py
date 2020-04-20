import warnings

import ipywidgets as widgets
import seaborn as sns
from SALib.analyze import sobol
from SALib.sample import saltelli
from ipywidgets import interact
from matplotlib import pyplot as plt

from .base_utils import _method_unit, _eprint
from .lca import *
from .lca import _expanded_names_to_names
from .params import _variable_params, _fixed_params, _param_registry


def _heatmap(df, title, vmax, ints=False):
    ''' Produce heatmap of a dataframe'''
    fig, ax = plt.subplots(figsize=(17, 17))
    sns.heatmap(df.transpose(), cmap="gist_heat_r", vmax=vmax, annot=True, fmt='.0f' if ints else '.2f', square=True)
    plt.title(title, fontsize=20)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)


def oat_matrix(model, impacts, n=10):
    '''Generates a heatmap of the incertitude of the model, varying input parameters one a a time.'''

    # Compile model into lambda functions for fast LCA
    lambdas, required_params = preMultiLCAAlgebric(model, impacts)

    required_param_names = _expanded_names_to_names(required_params)
    required_params = {key: _param_registry()[key] for key in required_param_names}
    var_params = _variable_params(required_param_names)


    change = np.zeros((len(var_params), len(impacts)))

    for iparam, param in enumerate(var_params.values()):
        params = {param.name: param.default for param in required_params.values()}

        # Compute range of values for given param
        params[param.name] = param.range(n)

        # Compute LCA
        df = postMultiLCAAlgebric(impacts, lambdas, **params)

        # Compute change
        change[iparam] = (df.max() - df.min()) / df.median() * 100

    # Build final heatmap
    change = pd.DataFrame(change, index=var_params.keys(), columns=[imp[2] for imp in impacts])
    _heatmap(change.transpose(), 'Change of impacts per variability of the input parameters (%)', 100, ints=True)


def _display_tabs(titlesAndContentF):
    '''Generate tabs'''
    tabs = []
    titles = []
    for title, content_f in titlesAndContentF:
        titles.append(title)

        tab = widgets.Output()
        with tab:
            content_f()
        tabs.append(tab)

    res = widgets.Tab(children=tabs)
    for i, title in enumerate(titles):
        res.set_title(i, title)
    display(res)


def oat_dasboard(modelOrLambdas, impacts, varying_param: ParamDef, n=10, all_param_names=None):
    '''
    Analyse the evolution of impacts for a single parameter. The other parameters are set to their default values.

    Parameters
    ----------
    model : activity, or lambdas as precomputed by preMultiLCAAlgebric, for faster computation
    impacts : set of methods
    param: parameter to analyse
    n: number of samples of the parameter
    '''

    if all_param_names == None:
        all_param_names = _param_registry().keys()

    params = {name: _param_registry()[name].default for name in all_param_names}

    # Compute range of values for given param
    params[varying_param.name] = varying_param.range(n)

    # print("Params: ", params)

    if isinstance(modelOrLambdas, Activity):
        df = multiLCAAlgebric(modelOrLambdas, impacts, **params)
    else:
        df = postMultiLCAAlgebric(impacts, modelOrLambdas, **params)

    # Add X values in the table
    pname = varying_param.name
    if varying_param.unit:
        pname = '%s [%s]' % (pname, varying_param.unit)
    df.insert(0, pname, varying_param.range(n))
    df = df.set_index(pname)

    def table():
        displayWithExportButton(df)

    def graph():

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            nb_rows = len(impacts) // 3 + 1

            fig, axes = plt.subplots(figsize=(15, 15))

            axes = df.plot(
                ax=axes, sharex=True, subplots=True,
                layout=(nb_rows, 3),
                # legend=None,
                kind='line' if varying_param.type == ParamType.FLOAT else 'bar')

            axes = axes.flatten()

            for ax, impact in zip(axes, impacts):
                ax.set_ylim(ymin=0)
                ax.set_ylabel(_method_unit(impact))

            plt.show(fig)

    def change():

        ch = (df.max() - df.min()) / df.median() * 100
        fig, ax = plt.subplots(figsize=(9, 6))
        plt.title('Relative change for %s' % df.index.name)
        ch.plot(kind='barh', rot=30)
        ax.set_xlabel('Relative change of the median value (%)')
        plt.tight_layout()
        plt.show(fig)

    _display_tabs([
        ("Graphs", graph),
        ("Data", table),
        ("Variation", change)
    ])


def oat_dashboard_interact(model, methods):
    '''Interative dashboard, with a dropdown for selecting parameter'''

    lambdas, required_params = preMultiLCAAlgebric(model, methods)

    required_params = _expanded_names_to_names(required_params)

    def process_func(param):
        oat_dasboard(lambdas, methods, _param_registry()[param], all_param_names=required_params)

    paramlist = list(_variable_params(required_params).keys())
    interact(process_func, param=paramlist)


def _stochastics(modelOrLambdas, methods, n=1000, var_params=None):
    ''' Compute stochastic impacts for later analysis of incertitude '''

    if var_params is None :
        var_params = _variable_params().values()

    # Extract variable names
    var_param_names = list([param.name for param in var_params])
    problem = {
        'num_vars': len(var_param_names),
        'names': var_param_names,
        'bounds': [[0, 1]] * len(var_param_names)
    }

    print("Generating samples ...")
    X = saltelli.sample(problem, n, calc_second_order=True)

    # Map normalized 0-1 random values into real values
    print("Transforming samples ...")
    params = dict()
    for i, param_name in enumerate(var_param_names):
        param = _param_registry()[param_name]
        params[param_name] = param.rand(X[:, i]).tolist()

    # Add static parameters
    for param in _param_registry().values():
        if param.name not in var_param_names :
            params[param.name] = param.default

    print("Processing LCA ...")
    if isinstance(modelOrLambdas, Activity):
        Y = multiLCAAlgebric(modelOrLambdas, methods, **params)
    else:
        Y = postMultiLCAAlgebric(methods, modelOrLambdas, **params)

    return problem, X, Y


def _sobols(methods, problem, Y):
    ''' Computes sobols indices'''
    s1 = np.zeros((len(problem['names']), len(methods)))
    s2 = np.zeros((len(problem['names']), len(problem['names']), len(methods)))
    st = np.zeros((len(problem['names']), len(methods)))

    def process(args) :
        imethod, method = args
        y = Y[Y.columns[imethod]]
        res = sobol.analyze(problem, y.to_numpy(), calc_second_order=True)
        return imethod, res

    with concurrent.futures.ThreadPoolExecutor() as exec :
        exec.map(process, enumerate(methods))

        for imethod, res in exec.map(process, enumerate(methods)):
            try:
                s1[:, imethod] = res["S1"]
                s2_ = np.nan_to_num(res["S2"])
                s2[:, :, imethod] = s2_ + np.transpose(s2_)
                st[:, imethod] = res["ST"]

            except Exception as e:
                _eprint("Sobol failed on %s" % imethod[2], e)

    return (s1, s2, st)





def _incer_stochastic_matrix(methods, param_names, Y, st):
    ''' Internal method computing matrix of parameter importance '''

    def draw(mode):

        if mode == 'sobol':
            data = st
        else:
            # If percent, express result as percentage of standard deviation / mean
            data = np.zeros((len(param_names), len(methods)))
            for i, method in enumerate(methods):
                # Total variance
                var = np.var(Y[Y.columns[i]])
                mean = np.mean(Y[Y.columns[i]])
                if mean != 0:
                    data[:, i] = np.sqrt((st[:, i] * var)) / mean * 100

        df = pd.DataFrame(data, index=param_names, columns=[method_name(method) for method in methods])
        _heatmap(
            df.transpose(),
            title="Relative deviation of impacts (%)" if mode == 'percent' else "Sobol indices (part of variability)",
            vmax=100 if mode == 'percent' else 1,
            ints=mode == 'percent')

    interact(draw, mode=[('Raw sobol indices (ST)', 'sobol'), ('Deviation (ST) / mean', 'percent')])


def incer_stochastic_matrix(modelOrLambdas, methods, n=1000, var_params=None):
    '''
    Method computing matrix of parameter importance

    parameters
    ----------
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    '''

    problem, X, Y = _stochastics(modelOrLambdas, methods, n, var_params)

    print("Processing Sobol indices ...")
    s1, s2, st = _sobols(methods, problem, Y)

    _incer_stochastic_matrix(methods, problem['names'], Y, st)


def _incer_stochastic_violin(methods, Y):
    ''' Internal method for computing violin graph of impacts '''

    nb_rows = math.ceil(len(methods) / 3)
    fig, axes = plt.subplots(nb_rows, 3, figsize=(15, 15), sharex=True)



    for imethod, method, ax in zip(range(len(methods)), methods, axes.flatten()):

        data = Y[Y.columns[imethod]]
        median = np.median(data)
        std = np.std(data)
        mean = np.mean(data)

        #ax.hist(Y[Y.columns[imethod]])
        ax.violinplot(data, showmedians=True)
        ax.title.set_text(method_name(method))
        ax.set_ylim(ymin=0)
        ax.set_ylabel(_method_unit(method))

        # Add text
        textstr = '\n'.join((
            r'$\mu=%.2f$' % (mean,),
            r'$\mathrm{median}=%.2f$' % (median,),
            r'$\sigma=%.2f$' % (std,)))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)



    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.show(fig)


def incer_stochastic_violin(modelOrLambdas, methods, n=1000, var_params=None):
    '''
    Method for computing violin graph of impacts

    parameters
    ----------
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    '''

    problem, X, Y = _stochastics(modelOrLambdas, methods, n, var_params)

    _incer_stochastic_violin(methods, Y)

percentiles = [10, 90, 25, 50, 75]
def _incer_stochastic_variations(methods, param_names, Y, sob1):
    ''' Method for computing violin graph of impacts '''
    method_names = [method_name(method) for method in methods]

    std = np.std(Y)
    mean = np.mean(Y)

    fig = plt.figure(num=None, figsize=(12, 6), dpi=80, facecolor='w', edgecolor='k')
    ax = plt.gca()
    tab20b = plt.get_cmap('tab20b')
    tab20c = plt.get_cmap('tab20c')
    ax.set_prop_cycle('color', [tab20b(k) if k < 1 else tab20c(k - 1) for k in np.linspace(0, 2, 40)])

    relative_variance_pct = std * std / (mean * mean) * 100
    totplt = plt.bar(np.arange(len(method_names)), relative_variance_pct, 0.8)

    sum = np.zeros(len(methods))

    plots = [totplt[0]]

    for i_param, param_name in enumerate(param_names):

        s1 = sob1[i_param, :]

        curr_bar = s1 * relative_variance_pct
        curr_plt = plt.bar(np.arange(len(method_names)), curr_bar, 0.8, bottom=sum)
        sum += curr_bar
        plots.append(curr_plt[0])

    plt.legend(plots, ['Higher order'] + param_names)
    plt.xticks(np.arange(len(method_names)), method_names, rotation=90)
    plt.title("variance / mean² (%)")
    plt.show(fig)

def _incer_stochastic_data(methods, param_names, Y, sob1, sobt):

    '''Show full stochastic output with sobol indices'''
    data = np.zeros((len(param_names) * 2 + len(percentiles) +2, len(methods)))
    data[0, :] = np.mean(Y)
    data[1, :] = np.std(Y)

    for i, percentile in enumerate(percentiles) :
        data[2 + i, :] = np.percentile(Y, percentile, axis=0)

    for i_param, param_name in enumerate(param_names):
        s1 = sob1[i_param, :]
        data[i_param + 2 + len(percentiles), :] = s1
        data[i_param + 2 + len(percentiles) + len(param_names), :] = sobt[i_param, :]

    rows = ["mean", "std"] + \
           ["p%d" % p for p in percentiles] + \
           ["Sobol 1(%s)" % param for param in param_names] + \
           ["Sobol T(%s)" % param for param in param_names]

    df = pd.DataFrame(data, index=rows, columns=[method_name(method) for method in methods])
    displayWithExportButton(df)

def incer_stochastic_dashboard(model, methods, n=1000, var_params=None):
    ''' Generates a dashboard with several statistics : matrix of parameter incertitude, violin diagrams, ...

    parameters
    ----------
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    '''

    problem, X, Y = _stochastics(model, methods, n, var_params)
    param_names = problem['names']

    print("Processing Sobol indices ...")
    s1, s2, st = _sobols(methods, problem, Y)

    def violin():
        _incer_stochastic_violin(methods, Y)

    def variation():
        _incer_stochastic_variations(methods, param_names, Y, s1)

    def matrix():
        _incer_stochastic_matrix(methods, problem['names'], Y, st)

    def data():
        _incer_stochastic_data(methods, problem['names'], Y, s1, st)

    _display_tabs([
        ("Violin graphs", violin),
        ("Impact variations", variation),
        ("Sobol matrix", matrix),
        ("Data", data)
    ])


def sobol_simplify_model(model, methods, min_ratio=0.8, n=1000, var_params=None) :
    '''
    Computes Sobol indices and selects main parameters for explaining sensibility of at least 'min_ratio',
    Then generates simplified models for thoses parameters.

    parameters
    ----------
    min_ratio: [0, 1] minimum amount of variation to explain
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    '''

    problem, X, Y = _stochastics(model, methods, n, var_params)
    param_names = problem['names']

    print("Processing Sobol indices ...")
    s1, s2, st = _sobols(methods, problem, Y)

    for imethod, method in enumerate(methods) :

        print("> Method : ", method_name(method))

        s1s2 = np.sum(s1[:, imethod]) + np.sum(s2[:, :, imethod]) / 2
        print("S1 + S2", s1s2)

        sum = 0
        sorted_param_indices = list(range(0, len(param_names)))
        sorted_param_indices = sorted(sorted_param_indices, key=lambda i : s1[i, imethod], reverse=True)
        selected_params = []
        for iparam, param in enumerate(sorted_param_indices) :

            selected_params.append(param_names[param])

            # S1
            sum += s1[param, imethod]

            # S2
            for iparam2 in range(0, iparam) :
                param2 = sorted_param_indices[iparam2]
                sum += s2[param, param2, imethod]

            if sum > min_ratio :
                break
        print("Selected params : ", selected_params, "explains: ", sum)

        fixedParams = [param for param in _param_registry().values() if param.name not in selected_params]

        # Generate simplified model
        simplified_models = simplifiedModel(
            model,
            [method],
            extraFixedParams=fixedParams)
        display(simplified_models[0])