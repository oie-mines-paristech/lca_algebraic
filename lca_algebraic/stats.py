import math
import random
import warnings
from time import time
from typing import Type, Dict, Tuple, List

import numpy as np
import seaborn as sns
from SALib.analyze import sobol
from SALib.sample import saltelli, sobol_sequence
from ipywidgets import interact
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from numpy import piecewise
from sympy import Float, Number, Add, AtomicExpr, Mul, Expr, Function, Abs, Sum, Eq, Piecewise
from sympy.core import symbol
from sympy.core.operations import AssocOp

from .base_utils import _method_unit
from .lca import *
from .lca import _expanded_names_to_names, _filter_param_values, _replace_fixed_params, _modelToExpr, _preMultiLCAAlgebric, _postMultiLCAAlgebric
from .params import _variable_params, _param_registry, FixedParamMode, _param_name, NameType

PARALLEL=False

def _parallel_map(f, items) :
    if PARALLEL :
        with concurrent.futures.ThreadPoolExecutor() as exec:
            return exec.map(f, items)
    else :
        return map(f, items)


def _heatmap(df, title, vmax, ints=False):
    ''' Produce heatmap of a dataframe'''
    fig, ax = plt.subplots(figsize=(17, 17))
    sns.heatmap(df.transpose(), cmap="gist_heat_r", vmax=vmax, annot=True, fmt='.0f' if ints else '.2f', square=True)
    plt.title(title, fontsize=20)
    plt.yticks(rotation=0)
    ax.tick_params(axis="x", labelsize=18)
    ax.tick_params(axis="y", labelsize=18)


def _extract_var_params(lambdas):
    required_param_names = set()
    for lamb in lambdas :
        required_param_names.update(_expanded_names_to_names(lamb.expanded_params))
    var_params = _variable_params(required_param_names)
    return sorted(var_params.values(), key=lambda p: (p.group if p.group else "", p.name))


@with_db_context(arg="model")
def oat_matrix(model, impacts, functional_unit=1, n=10, title='Impact variability (% of mean)', name_type=NameType.LABEL):
    '''Generates a heatmap of the incertitude of the model, varying input parameters one a a time.'''

    # Compile model into lambda functions for fast LCA
    lambdas = _preMultiLCAAlgebric(model, impacts, alpha=1/functional_unit)

    # Sort params by category
    sorted_params = _extract_var_params(lambdas)

    change = np.zeros((len(sorted_params), len(impacts)))

    for iparam, param in enumerate(sorted_params):
        params = {param.name: param.default for param in sorted_params}

        # Compute range of values for given param
        params[param.name] = param.range(n)

        # Compute LCA
        df = _postMultiLCAAlgebric(impacts, lambdas, **params)

        # Compute change
        change[iparam] = (df.max() - df.min()) / df.median() * 100



    # Build final heatmap
    change = pd.DataFrame(change,
                index=[_param_name(param, name_type) for param in sorted_params],
                columns=[method_name(imp) for imp in impacts])
    _heatmap(change.transpose(), title, 100, ints=True)


def _display_tabs(titlesAndContentF):
    """Generate tabs"""
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


def oat_dasboard(modelOrLambdas, impacts, varying_param: ParamDef, n=10, all_param_names=None,
                 figsize=(15, 15), figspace=(0.5, 0.5), sharex=True, cols=3, func_unit="kWh"):
    '''
    Analyse the evolution of impacts for a single parameter. The other parameters are set to their default values.
    The result heatmap shows percentage of variation relative to median value.

    Parameters
    ----------
    model : activity, or lambdas as precomputed by preMultiLCAAlgebric, for faster computation
    impacts : set of methods
    param: parameter to analyse
    n: number of samples of the parameter
    figsize: Size of figure fro graphs : (15, 15 by default)
    figspace: Space between figures for graphs : (0.5, 0.5) by default
    sharex: Shared X axes ? True by default
    '''

    if all_param_names == None:
        all_param_names = _param_registry().keys()

    params = {name: _param_registry()[name].default for name in all_param_names}

    # Compute range of values for given param
    params[varying_param.name] = varying_param.range(n)

    # print("Params: ", params)

    if isinstance(modelOrLambdas, Activity):
        with DbContext(modelOrLambdas) :
            df = compute_impacts(modelOrLambdas, impacts, **params)
    else:
        df = _postMultiLCAAlgebric(impacts, modelOrLambdas, **params)

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

            nb_rows = int(math.ceil(len(impacts) / cols))

            fig, axes = plt.subplots(figsize=figsize)
            plt.subplots_adjust(None, None, None, None, figspace[0], figspace[1])

            axes = df.plot(
                ax=axes, sharex=sharex, subplots=True,
                layout=(nb_rows, cols),
                # legend=None,
                kind='line' if varying_param.type == ParamType.FLOAT else 'bar')

            axes = axes.flatten()

            for ax, impact in zip(axes, impacts):
                ax.set_ylim(ymin=0)
                unit = _method_unit(impact) + " / " + func_unit
                ax.set_ylabel(unit, fontsize=15)
                ax.set_xlabel(pname, fontsize=15)

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

@with_db_context(arg="model")
def oat_dashboard_interact(model, methods, functional_unit=1, **kwparams):
    '''Interactive dashboard, with a dropdown for selecting parameter

    Parameters
    ----------
    figsize: Size of figure fro graphs : (15, 15 by default)
    figspace: Space between figures for graphs : (0.5, 0.5) by default
    sharex: Shared X axes ? True by default
    '''

    lambdas = _preMultiLCAAlgebric(model, methods, alpha=1/functional_unit)

    def process_func(param):
        with DbContext(model):
            oat_dasboard(lambdas, methods, _param_registry()[param], **kwparams)

    param_list = _expanded_names_to_names(lambdas[0].expanded_params)
    param_list = list(_variable_params(param_list).keys())

    interact(process_func, param=param_list)


class StochasticMethod :
    SALTELLI = "saltelli"
    RAND = "rand"
    SOBOL="sobol"


def _stochastics(
        modelOrLambdas, methods, n=1000,
        var_params=None,
        sample_method=StochasticMethod.SALTELLI,
        functional_unit=1,
        **extra_fixed_params):

    params, problem = _generate_random_params(n, sample_method, var_params)

    # Fix other params
    if extra_fixed_params :
        params.update(extra_fixed_params)

    Y = _compute_stochastics(modelOrLambdas, methods, params=params, functional_unit=functional_unit)

    return problem, params, Y


def _compute_stochastics(modelOrLambdas, methods, functional_unit=1, params=dict()):
    if isinstance(modelOrLambdas, Activity):
        Y = compute_impacts(modelOrLambdas, methods, functional_unit=functional_unit, **params)
    else:
        Y = _postMultiLCAAlgebric(methods, modelOrLambdas, **params)
    return Y

def _generate_random_params(n, sample_method=StochasticMethod.SALTELLI, var_params=None, seed=None):

    ''' Compute stochastic impacts for later analysis of incertitude '''
    if var_params is None:
        var_params = _variable_params().values()

    if seed is None :
        seed = int(time() * 1000)

    random.seed(seed)

    # Extract variable names
    var_param_names = list([param if isinstance(param, str) else param.name for param in var_params])
    problem = {
        'num_vars': len(var_param_names),
        'names': var_param_names,
        'bounds': [[0, 1]] * len(var_param_names)}
    print("Generating samples ...")
    if sample_method == StochasticMethod.SALTELLI:
        X = saltelli.sample(problem, n, calc_second_order=True)
    elif sample_method == StochasticMethod.RAND:
        X = np.random.rand(n, len(var_param_names))
    elif sample_method == StochasticMethod.SOBOL:
        X = sobol_sequence.sample(n * (len(var_param_names) * 2 + 2), len(var_param_names))
    # elif sample_method == StochasticMethod.LATIN :
    #    X = latin.sample(problem, n)
    else:
        raise Exception("Unkown rand method " + sample_method)
    # Map normalized 0-1 random values into real values
    print("Transforming samples ...")
    params = dict()
    for i, param_name in enumerate(var_param_names):
        param = _param_registry()[param_name]
        params[param_name] = param.rand(X[:, i]).tolist()

    # Add fixed parameters
    for param in _param_registry().values():
        if param.name not in var_param_names:
            params[param.name] = param.default

    return params, problem


class SobolResults :

    def __init__(self, s1, s2, st, s1_conf=None, s2_conf=None, st_conf=None):
        self.s1 = s1
        self.s2 = s2
        self.st = st
        self.s1_conf = s1_conf
        self.s2_conf = s2_conf
        self.st_conf = st_conf


def _sobols(methods, problem, Y) -> SobolResults :
    ''' Computes sobols indices'''
    s1 = np.zeros((len(problem['names']), len(methods)))
    s1_conf = np.zeros((len(problem['names']), len(methods)))
    s2 = np.zeros((len(problem['names']), len(problem['names']), len(methods)))
    s2_conf = np.zeros((len(problem['names']), len(problem['names']), len(methods)))
    st = np.zeros((len(problem['names']), len(methods)))
    st_conf = np.zeros((len(problem['names']), len(methods)))

    def process(args) :
        imethod, method = args

        print("Processing sobol for " + str(method))
        y = Y[Y.columns[imethod]]
        res = sobol.analyze(problem, y.to_numpy(), calc_second_order=True)
        return imethod, res

    for imethod, res in _parallel_map(process, enumerate(methods)):
        try:
            s1[:, imethod] = res["S1"]
            s1_conf[:, imethod] = res["S1_conf"]
            s2_ = np.nan_to_num(res["S2"])
            s2_conf_ = np.nan_to_num(res["S2_conf"])
            s2[:, :, imethod] = s2_ + np.transpose(s2_)
            s2_conf[:, :, imethod] = s2_conf_ + np.transpose(s2_conf_)
            st[:, imethod] = res["ST"]
            st_conf[:, imethod] = res["ST_conf"]

        except Exception as e:
            error("Sobol failed on %s" % imethod[2], e)

    return SobolResults(s1, s2, st, s1_conf, s2_conf, st_conf)





def _incer_stochastic_matrix(methods, param_names, Y, sob, name_type=NameType.LABEL):
    ''' Internal method computing matrix of parameter importance '''



    def draw(indice, mode):
        '''
        Mode comes as ('s1' || 'st', 'raw' || 'percent')
        '''

        sx = sob.s1 if indice == "s1" else sob.st
        if mode == 'raw':
            data = sx
        else:
            # If percent, express result as percentage of standard deviation / mean
            data = np.zeros((len(param_names), len(methods)))
            for i, method in enumerate(methods):
                # Total variance
                var = np.var(Y[Y.columns[i]])
                mean = np.mean(Y[Y.columns[i]])
                if mean != 0:
                    data[:, i] = np.sqrt((sx[:, i] * var)) / mean * 100

        param_labels = [_param_name(_param_registry()[name], name_type) for name in param_names]
        df = pd.DataFrame(data, index=param_labels, columns=[method_name(method) for method in methods])
        _heatmap(
            df.transpose(),
            title="Relative deviation of impacts (%)" if mode == 'percent' else "Sobol indices (part of variability)",
            vmax=100 if mode == 'percent' else 1,
            ints=mode == 'percent')

    interact(draw,
             indice=['s1', 'st'],
             mode=[
                 ('Raw indices', 'raw'),
                 ('Relative to mean (%)', 'percent')]
             )

@with_db_context(arg="model")
def incer_stochastic_matrix(model, methods, functional_unit=1, n=1000, name_type=NameType.LABEL):
    '''
    Method computing matrix of parameter importance

    parameters
    ----------
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    '''

    lambdas = _preMultiLCAAlgebric(model, methods, alpha=1/functional_unit)
    var_params = _extract_var_params(lambdas)

    problem, _, Y = _stochastics(lambdas, methods, n, var_params)

    print("Processing Sobol indices ...")
    sob = _sobols(methods, problem, Y)

    _incer_stochastic_matrix(methods, problem['names'], Y, sob, name_type=name_type)


def _incer_stochastic_violin(methods, Y, figsize=(15, 15), figspace=(0.5, 0.5), sharex=True, nb_cols=3):
    ''' Internal method for computing violin graph of impacts
    Parameters
    ----------
    methods: list of impact methods
    Y : output
    figsize: Size of figure for graphs : (15, 15 by default)
    figspace: Space between figures for graphs : (0.5, 0.5) by default
    sharex: Shared X axes ? True by default
    nb_cols: Number of colums. 3 by default
    '''

    nb_rows = math.ceil(len(methods) / nb_cols)
    fig, axes = plt.subplots(nb_rows, nb_cols, figsize=figsize, sharex=sharex)
    plt.subplots_adjust(None, None, None, None, figspace[0], figspace[1])

    for imethod, method, ax in zip(range(len(methods)), methods, axes.flatten()):

        data = Y[Y.columns[imethod]]
        median = np.median(data)
        std = np.std(data)
        mean = np.mean(data)

        ax.violinplot(data, showmedians=True)
        ax.title.set_text(method_name(method))
        ax.set_ylim(ymin=0)
        ax.set_ylabel(_method_unit(method))
        ax.set_xticklabels([])

        # Add text
        textstr = '\n'.join((
            r'$\mu=%.3g$' % (mean,),
            r'$\mathrm{median}=%.3g$' % (median,),
            r'$\sigma=%.3g$' % (std,)))
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
                verticalalignment='top', bbox=props)



    plt.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    plt.show(fig)


@with_db_context(arg="modelOrLambdas")
def incer_stochastic_violin(modelOrLambdas, methods, functional_unit=1, n=1000, var_params=None, **kwparams):
    '''
    Method for computing violin graph of impacts

    parameters
    ----------
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    '''
    _, _, Y = _stochastics(
        modelOrLambdas, methods, n=n,
        var_params=var_params,
        functional_unit=functional_unit)

    _incer_stochastic_violin(methods, Y, **kwparams)

_percentiles = [10, 90, 25, 50, 75]
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

    plt.legend(plots, ['Higher order'] + param_names, loc=(1, 0))
    plt.xticks(np.arange(len(method_names)), method_names, rotation=90)
    plt.title("variance / mean² (%)")
    plt.show(fig)

def _incer_stochastic_data(methods, param_names, Y, sob1, sobt):

    '''Show full stochastic output with sobol indices'''
    data = np.zeros((len(param_names) * 2 + len(_percentiles) +2, len(methods)))
    data[0, :] = np.mean(Y)
    data[1, :] = np.std(Y)

    for i, percentile in enumerate(_percentiles) :
        data[2 + i, :] = np.percentile(Y, percentile, axis=0)

    for i_param, param_name in enumerate(param_names):
        s1 = sob1[i_param, :]
        data[i_param + 2 + len(_percentiles), :] = s1
        data[i_param + 2 + len(_percentiles) + len(param_names), :] = sobt[i_param, :]

    rows = ["mean", "std"] + \
           ["p%d" % p for p in _percentiles] + \
           ["Sobol 1(%s)" % param for param in param_names] + \
           ["Sobol T(%s)" % param for param in param_names]

    df = pd.DataFrame(data, index=rows, columns=[method_name(method) for method in methods])
    displayWithExportButton(df)

@with_db_context(arg="model")
def incer_stochastic_dashboard(model, methods, n=1000, var_params=None, functional_unit=1, **kwparams):
    ''' Generates a dashboard with several statistics : matrix of parameter incertitude, violin diagrams, ...

    parameters
    ----------
    var_params: Optional list of parameters to vary.
    By default use all the parameters with distribution not FIXED
    figsize: Size of figure for violin plots : (15, 15) by default
    figspace: Space between violin graphs (0.5, 0.5) by default
    sharex: Share X axe for violin graph : True by default
    '''

    problem, _, Y = _stochastics(
        model, methods, n,
        var_params=var_params, functional_unit=functional_unit)

    param_names = problem['names']

    print("Processing Sobol indices ...")
    sob = _sobols(methods, problem, Y)

    def violin():
        _incer_stochastic_violin(methods, Y, **kwparams)

    def variation():
        _incer_stochastic_variations(methods, param_names, Y, sob.s1)

    def matrix():
        _incer_stochastic_matrix(methods, problem['names'], Y, sob)

    def data():
        _incer_stochastic_data(methods, problem['names'], Y, sob.s1, sob.st)

    _display_tabs([
        ("Violin graphs", violin),
        ("Impact variations", variation),
        ("Sobol matrix", matrix),
        ("Data", data)
    ])


def _round_expr(expr, num_digits):
    ''' Round all number in sympy expression with n digits'''

    return expr.xreplace({n : Float(n, num_digits) if isinstance(n, Float) else n for n in expr.atoms(Number)})

def _snake2camel(val) :
    return ''.join(word.title() for word in val.split('_'))


def _enum_to_piecewize(exp) :

    def checkEnumSymbol(term) :
        if not isinstance(term, Symbol) :
            return (None, None)
        if not "_" in term.name :
            return (None, None)
        enum_name, enum_val = term.name.rsplit("_", 1)

        if not enum_name in _param_registry() :
            return (None, None)
        return enum_name, enum_val

    def checkEnumProduct(term) :

        """If term is enumVal * X return (param, value, X) else return null"""

        if not isinstance(term, Mul) :
            return (None, None, None)

        if len(term.args) != 2 :
            return (None, None, None)

        a, b = term.args

        name, value = checkEnumSymbol(a)
        if name is not None :
            return (name, value, b)

        name, value = checkEnumSymbol(b)
        if name is not None:
            return (name, value, a)

        return (None, None, None)


    def _replace_enums(expr) :
        # Dict of enum_name -> { enum_value -> ratio }
        enums = defaultdict(lambda : dict())
        res_terms = []

        for term in expr.args :

            name, value, ratio = checkEnumProduct(term)
            if name is not None :
                # This is a enum value !
                enums[name][value] = ratio
            else :
                # Other term
                res_terms.append(term)

        if len(enums) == 0 :
            # No change
            return expr

        for enum_name, ratio_dict in enums.items() :
            choices = [(ratio, Eq(symbols(enum_name), symbols(enum_value))) for enum_value, ratio in ratio_dict.items()]
            if len(choices) < len(_param_registry()[enum_name].values) :
                # Not all choices covered ? => Add default
                choices.append((0, True))

            res_terms.append(Piecewise(*choices))

        return Add(*res_terms)

    return exp.replace(
        lambda x: isinstance(x, Sum) or isinstance(x, Add),
        _replace_enums)

def prettify(exp) :
     """
     Prettify expression for publication :
     > change snake_symbols to SnakeSymbols (avoiding lowerscript in Latex)"""

     res = _enum_to_piecewize(exp)

     res = res.replace(
         lambda x: isinstance(x, Symbol),
         lambda x: Symbol(_snake2camel(str(x))))

     # Replace absolute values for positive parameters
     res, nb_match = _replace_abs(res)
     if nb_match > 0:
         # It changed => simplify again
         res = simplify(res)

     return res

def _replace_abs(exp) :
    """Replace |X| by X if X is float param """

    nb_match = 0

    def replaceAbs(absExp : Abs) :
        nonlocal nb_match
        if len(absExp.args) != 1 :
            return absExp
        arg = absExp.args[0]

        if not isinstance(arg, Symbol) :
            return absExp

        params = _param_registry()
        if not arg.name in params :
            return absExp

        param = params[arg.name]
        if param.type == ParamType.FLOAT and param.min >= 0 :
            nb_match += 1
            return arg
        else :
            return absExp

    res = exp.replace(
        lambda x: isinstance(x, Abs),
        lambda x: replaceAbs(x))

    return res, nb_match

@with_db_context(arg="model")
def sobol_simplify_model(
    model, methods,
    min_ratio=0.8, n=2000, var_params=None,
    fixed_mode = FixedParamMode.MEDIAN,
    num_digits=3,
    simple_sums=True,
    functional_unit=1,
    simple_products=True) -> List[LambdaWithParamNames]:

    '''
    Computes Sobol indices and selects main parameters for explaining sensibility of at least 'min_ratio',
    Then generates simplified models for those parameters.

    parameters
    ----------
    min_ratio: [0, 1] minimum amount of first order variation (sum of S1) to explain
    var_params: Optional list of parameters to vary.
    fixed_mode : What to replace minor parameters with : MEDIAN by default
    sob: [optional] Pre-computed sobol indices
    simplify_sums: If true (default) remove terms in sums that are lower than 1%

    returns
    _______

    List of LambdaWithParamNames, one per impact : with wraps the simplified expression together with the
    list of required parameters and a fast complied lambda function for fast evaluation.
    '''

    # Default var param names
    if var_params is None :
        var_params = _variable_params().values()

    var_param_names = list([param.name for param in var_params])

    problem, params, Y = _stochastics(
        model, methods, n,
        var_params=var_params,
        functional_unit=functional_unit)

    sob = _sobols(methods, problem, Y)

    s1, s2 = sob.s1, sob.s2

    res = []

    # Generate simplified model
    exprs = _modelToExpr(
        model, methods,
        alpha=1/functional_unit)

    for imethod, method in enumerate(methods) :

        print("> Method : ", method_name(method))

        s1_sum = np.sum(s1[:, imethod])
        s2_sum = np.sum(s2[:, :, imethod]) / 2
        print('S1: ', s1_sum)
        print('S2: ', s2_sum)
        print('ST: ', np.sum(sob.st[:, imethod]))

        sum = 0

        sorted_param_indices = list(range(0, len(var_param_names)))
        sorted_param_indices = sorted(sorted_param_indices, key=lambda i : s1[i, imethod], reverse=True)
        selected_params = []
        sobols = dict()

        for iparam, param in enumerate(sorted_param_indices) :


            # S1
            sum += s1[param, imethod]

            param_name = var_param_names[param]
            selected_params.append(param_name)
            sobols[param_name] = s1[param, imethod]

            # S2
            #for iparam2 in range(0, iparam) :
            #    param2 = sorted_param_indices[iparam2]
            #    sum += s2[param, param2, imethod]

            if sum > min_ratio :
                break
        print("Selected params : ", selected_params, "explains: ", sum)

        expr = exprs[imethod]

        # Replace non selected params by their value
        fixed_params = [param for param in _param_registry().values() if param.name not in selected_params]
        expr = _replace_fixed_params(expr, fixed_params, fixed_mode=fixed_mode)

        # Sympy simplification
        expr = simplify(expr)

        # Round numerical values to 3 digits
        expr = _round_expr(expr, num_digits)

        # Lambdify the expression
        lambd = LambdaWithParamNames(expr, params=selected_params, sobols=sobols)

        # Compute list of parameter values (monte carlo)
        completed_params = lambd.complete_params(params)
        expanded_params = _filter_param_values(completed_params, lambd.expanded_params)

        # Extra step of simplification : simplify sums with neligeable terms
        if simple_sums :
            expr = _simplify_sums(expr, expanded_params)

        if simple_products:
            expr = _simplify_products(expr, expanded_params)


        expr = simplify(expr)

        display(prettify(expr))

        res.append(LambdaWithParamNames(expr, params=selected_params, sobols=sobols))

    return res

TERM_MIN_LEVEL = 0.01

def _rec_expression(exp, func) :
    """ Recurse trough an expression, transforming each term with the result of f(term) """
    def rec(term) :
        if issubclass(exp.func, AtomicExpr):
            return func(term)
        else :
            args = filter(lambda x: x is not None, list(func(arg) for arg in term.args))
            return term.func(*args)

def _simplify_sums(expr, param_values) :

    def replace_term(term, minv, maxv, max_max):
        abs_max = max(abs(minv), abs(maxv))
        if abs_max < (TERM_MIN_LEVEL * max_max):
            return None
        else:
            return term

    return _simplify_terms(expr, param_values, Add, replace_term)

def _simplify_products(expr, param_values) :

    def replace_term(term, minv, maxv, max_max) :

        # Close to 1 or -1 ?
        for factor in [-1, 1]:
            if abs(minv - factor) < TERM_MIN_LEVEL and abs(maxv - factor) < TERM_MIN_LEVEL:
                if factor == -1:
                    return -1
                else:
                    # * 1.0 : remove term
                    return None

        return term

    return _simplify_terms(expr, param_values, Mul, replace_term)

def _simplify_terms(expr, expanded_param_values, op:Type[AssocOp], replace) :

    """Generic simplification of sum or product"""

    # Determine max normalized value of this term, for all param values (monte carlo)
    min_max_cache : Dict[str, Tuple[float, float]] = dict()
    def min_max(term) :

        # In cache ?
        key = str(term)
        if key in min_max_cache :
            return min_max_cache[key]

        # Non varying ?
        if len(term.free_symbols) == 0:
            values = [term.evalf()]
        else:
            lambd_term = lambdify(expanded_param_values.keys(), term)
            values = lambd_term(**expanded_param_values)

        minv = np.min(values)
        maxv = np.max(values)
        min_max_cache[key] = (minv, maxv)
        return (minv, maxv)

    # Cleanup :keep only most relevant terms
    def cleanup(exp):

        if (not isinstance(exp, Expr)) or issubclass(exp.func, AtomicExpr):
            return exp

        # For Op, only select terms than are relevant
        if exp.func == op :

            # Compute max of max
            def abs_max(minv, maxv) :
                return max(abs(minv), abs(maxv))

            max_max = max([abs_max(*min_max(arg)) for arg in exp.args])

            # Only keep term above level
            args = [replace(arg, *min_max(arg), max_max) for arg in exp.args]
        else:
            args = exp.args

        args = [cleanup(arg) for arg in args if arg is not None]

        return exp.func(*args)

    return cleanup(expr)



def _hline(x1, x2, y, linewidth=1, linestyle='solid'):
    ymin, ymax = plt.ylim()
    xmin, xmax = plt.xlim()
    minx = (x1 - xmin) / (xmax - xmin)
    maxx = (x2 - xmin) / (xmax - xmin)
    plt.axhline(ymax * y, color='k', xmin=minx, xmax=maxx, linewidth=linewidth, linestyle=linestyle)


def _vline(x, ymin, ymax, linewidth=1, linestyle='solid'):
    plt.axvline(x, color='k', ymin=ymin, ymax=ymax, linewidth=linewidth, linestyle=linestyle)


def _graph(data, unit, title, ax, alpha=1, textboxtop=0.95, textboxright=0.95, color=None,
           limit_xrange=False, percentiles=[5,95], fontsize=12):


    if ax is not None:
        plt.sca(ax)
    else:
        ax = plt.gca()

    median = np.median(data)
    std = np.std(data)
    mean = np.mean(data)
    xmin= np.min(data)

    p9995 = np.percentile(data, 99.95)

    pvals = [np.percentile(data, perc) for perc in percentiles]

    variability = std / mean

    args = dict()
    if color:
        args['color'] = color

    if limit_xrange :
        plt.xlim(xmin, p9995)
    plt.hist(data, 200, alpha=alpha, **args)

    perc_strs = [r'$p%d=%.3g$' % (p,pval) for p, pval in zip(percentiles, pvals)]

    textstr = '\n'.join([
        r'$\mu=%.3g$' % (mean,),
        r'$\mathrm{median}=%.3g$' % (median,),
        r'$\sigma=%.3g$' % (std,),
        r'$\sigma/\mu=%.3g$' % (variability,)
    ] + perc_strs)

    props = dict(boxstyle='round', facecolor='wheat' if not color else color, alpha=0.5)
    ax.text(textboxright, textboxtop, textstr, transform=ax.transAxes, fontsize=fontsize,
            verticalalignment='top', ha='right', bbox=props)

    # Axes
    ax.set_xlabel(unit, dict(fontsize=14))
    ax.set_yticks([])
    ax.set_title(title, dict(fontsize=16))

    return dict(
        median=median,
        std=std,
        p=pvals,
        mean=mean,
        var=variability)

def distrib(*args, **kwargs) :
    """
       Show distributions together with statistical outcomes

       Synonym of #graphs()

       parameters
       ----------
       model: normalized model
       methods: List of impacts
       Y: output of processing. If None, monte carlo will be processed again
       nb_cols : number of colons to display graphs on
       invert : list of methods for which result should be inverted (1/X). None by default
       scales : Dict of method => scale, for multiplying results. To be used with unit overrides
       unit_overrides : Dict of method => string for overriding unit, in respect to custom scales
       height: Height of graph : 10 inches be default
       width : Width of graphs : 15 inches by default
       """
    return graphs(*args, **kwargs)


@with_db_context(arg="model")
def graphs(
        model, methods,
        functional_unit=1,
        Y=None, nb_cols=1, axes=None, title=None,
        invert=None,
        scales=None, # Dict of method => scale
        unit_overrides=None,
        height=10, width=15,
        func_unit="kWh",
        **kwargs):
    """
    Show distributions together with statistical outcomes

    parameters
    ----------
    model: normalized model
    methods: List of impacts
    Y: output of processing. If None, monte carlo will be processed again
    nb_cols : number of colons to display graphs on
    invert : list of methods for which result should be inverted (1/X). None by default
    scales : Dict of method => scale, for multiplying results. To be used with unit overrides
    unit_overrides : Dict of method => string for overriding unit, in respect to custom scales
    height: Height of graph : 10 inches be default
    width : Width of graphs : 15 inches by default
    """


    if Y is None:
        _, _, Y = _stochastics(model, methods, n=10000, functional_unit=functional_unit)

    if axes is None:
        nb_rows = math.ceil(len(methods) / nb_cols)
        fig, axes = plt.subplots(nb_rows, nb_cols, figsize=(width, height * nb_rows))

    if isinstance(axes, np.ndarray):
        axes = axes.flatten()
    else:
        axes = [axes]

    plt.subplots_adjust(hspace=0.4)


    res=dict()
    for i, method, ax in zip(range(len(methods)), methods, axes):

        data = Y[Y.columns[i]]

        if invert and method in invert:
            data = 1 / data

        if scales and method in scales:
            data = data * scales[method]

        graph_title = title if title else method_name(method)

        if unit_overrides and method in unit_overrides:
            unit = unit_overrides[method]
        else:
            unit = _method_unit(method)

        unit += " / " + func_unit

        stats = _graph(
            data, unit,
            graph_title,
            ax=ax,
            **kwargs)

        res[graph_title + (' [%s]' % unit)] = stats

    for i in range(0, -len(methods) % nb_cols):
        ax = axes.flatten()[-(i + 1)]
        ax.axis("off")

    return pd.DataFrame(res)


@with_db_context(arg="model")
def compare_simplified(
        model, methods, simpl_lambdas,
        functional_unit=1,
        scales=None,  # Dict of method => scale
        unit_overrides=None,
        nb_cols=2, height=10, width=15, textboxright=0.6, r2_height=0.65, func_unit="kWh",
        residuals=False,
        **kwargs):
    '''
    Compare distribution of simplified model with full model

    Parameters
    ----------
    model: Model
    residuals : If true, draw heat map of residuals, instead of distributions
    methods : Impact methods
    simpl_lambdas : Simplified lambdas, as returned by sobol_simplify_model(...)
    nb_cols: number of columns for displaying graphs
    percentiles: List of percentiles to compute [5, 95] by default
    '''

    # Raw model
    lambdas = _preMultiLCAAlgebric(
        model,
        methods,
        alpha=1/functional_unit)

    nb_rows = math.ceil(len(methods) / nb_cols)
    fig, axes = plt.subplots(nb_rows, nb_cols, figsize=(width, height * nb_rows))

    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3)

    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for i, lambd, simpl_lambd, method, ax in zip(range(len(methods)), lambdas, simpl_lambdas, methods, axes.flatten()):
        params, _ = _generate_random_params(100000, sample_method=StochasticMethod.RAND)

        # Run  Monte Carlo on full model
        Y1 = _compute_stochastics([lambd], [method], params)
        d1 = Y1[Y1.columns[0]]

        # Run monte carlo of simplified model
        Y2 = _compute_stochastics([simpl_lambd], [method], params)
        d2 = Y2[Y2.columns[0]]

        r_value = r_squared(Y1, Y2)

        title = method_name(method)

        if scales and method in scales:
            d1 = d1 * scales[method]
            d2 = d2 * scales[method]

        if unit_overrides and method in unit_overrides:
            unit = unit_overrides[method]
        else:
            unit = _method_unit(method)

        unit += " / " + func_unit

        if residuals :
            hb = ax.hexbin(d1, d2, gridsize=(200, 200), mincnt=1)
            cb = fig.colorbar(hb, ax=ax)

            xmin, xmax = ax.get_xlim()
            ymin, ymax = ax.get_ylim()
            vmax= max(ymax, xmax)
            line = Line2D([0, vmax], [0, vmax], color='grey', linestyle='dashed')
            ax.add_line(line)
            ax.set_xlim(0, vmax)
            ax.set_ylim(0, vmax)
            cb.set_label("Counts")
            ax.set_xlabel("Reference model [%s]" % (unit), dict(fontsize=14))
            ax.set_ylabel("Simplified model [%s]" % (unit), dict(fontsize=14))
            ax.set_title(title, dict(fontsize=16))

        else:
            _graph(d1, unit, title, ax=ax, alpha=0.6, color=colors[0], **kwargs)
            _graph(d2, unit, title, ax=ax, alpha=0.6, textboxright=textboxright, color=colors[1], **kwargs)

        ax.text(0.9, r2_height, "R² : %0.3g" % r_value, transform=ax.transAxes, fontsize=14,
                verticalalignment='top', ha='right')

    # Hide missing graphs
    for i in range(0, -len(methods) % nb_cols):
        ax = axes.flatten()[-(i + 1)]