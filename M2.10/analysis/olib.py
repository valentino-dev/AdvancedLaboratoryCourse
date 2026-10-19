#from pathvalidate import sanitize_filename
import os
#import pdfkit
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.differentiate import jacobian
from functools import partial

def get_sigmas(model, X, popt, pcov):
    j = np.array([jacobian(lambda p: partial(model, x)(*p), popt).df for x in X])
    sigma = np.sqrt(np.einsum('ij,jk,ik->i', j, pcov, j))
    return sigma

def calc_reduced_chi_sq(model, x, y, y_err, popt):
    residuals = y-model(x, *popt)
    chi_sq = ((residuals/y_err)**2).sum()
    dof = len(y) - len(popt)
    return chi_sq / dof

def fit_function(model, table, p0, plot=False, k=1, show_p0=False, label='', color='red', show_band=True):
    x, y, y_err = table.x, table.y, table.y_err
    if np.all(y_err):
        popt, pcov = curve_fit(model, x, y, p0=p0, sigma=y_err, absolute_sigma=True)
    else:
        popt, pcov = curve_fit(model, x, y, p0=p0)
    returns = [popt, pcov]
    red_chi_sq = calc_reduced_chi_sq(model, x, y, y_err, popt)
    print(f"$\chi^2/\\text{{ddof}}: \\num{{%.2f}}$"%(red_chi_sq))
    for p, p_err in zip(popt, np.sqrt(np.diag(pcov))):
        p, p_err = roundToError(p, p_err)
        print(f'\SI{{{p}+-{p_err}}}{{}}')

    if plot:
        X = np.linspace(min(x), max(x), 200)
        if show_p0:
            Y = model(X, *p0)
        else:
            Y = model(X, *popt)
            sigma = get_sigmas(model, X, popt, pcov)
            returns.append(sigma)
            if show_band:
                plt.fill_between(X/table.x_scaling, (Y-k*sigma)/table.y_scaling, (Y+k*sigma)/table.y_scaling, alpha=0.3, color=color)
        plt.plot(X/table.x_scaling, Y/table.y_scaling, color=color, label=label+f'Fit: $\\chi^2/\\text{{ddof}}={red_chi_sq:.2f}$')

    return returns


def roundToError(x, x_err, digits=2):
    if type(x) == np.ndarray:
        error_digits = np.zeros(x_err.shape, dtype=np.int32)
        error_digits[x_err!=0. * ~np.isinf(x_err)] = np.array(np.floor(np.log10(x_err[x_err!=0.*~np.isinf(x_err)])), dtype=np.int32)
        for i in range(x.shape[0]):
            x[i]=np.round(x[i], digits-error_digits[i]-1)
            x_err[i]=np.round(x_err[i], digits-error_digits[i]-1)

        return x, x_err
    elif x_err != 0 and not math.isinf(x_err):
        error_digits = int(np.floor(np.log10(x_err)))
        return np.round(x, digits-error_digits-1), np.round(x_err, digits-error_digits-1)
    else:
        return x, x_err


class model:
    def __init__(self) -> None:
        pass

    def fit(self, X, Y):
        pass

    def predict(self, x) -> np.ndarray:
        pass


class Table:

    def __init__(self, x, x_err, y, y_err, title=None, xlabel="x", ylabel="y", xunit='a.u.', yunit='a.u.'):
        self.x = x
        self.x_err = x_err
        self.y = y
        self.y_err = y_err
        
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.xunit = xunit
        self.yunit = yunit

        self.change_scale(*self.getScaling())
        

    def change_scale(self, x_scaling, y_scaling):
        self.x_scaling, self.y_scaling = x_scaling, y_scaling

        self.x_prefix, self.y_prefix = self.get_prefix(self.x_scaling), self.get_prefix(self.y_scaling)
        self.x_unit_scaled, self.y_unit_scaled = self.x_prefix+self.xunit, self.y_prefix+self.yunit

        self.x_scaled = self.x / self.x_scaling
        self.x_err_scaled = self.x_err / self.x_scaling

        self.y_scaled = self.y / self.y_scaling
        self.y_err_scaled = self.y_err / self.y_scaling


    def get_prefix(self, scaling):
        magnitude = int(np.log10(scaling))//3
        if magnitude == 0:
            return ''
        elif magnitude > 0 and magnitude < 6:
            return ['k', 'M', 'G', 'T', 'P'][magnitude-1]
        elif -magnitude > 0 and -magnitude < 6:
            return ['m', 'µ', 'n', 'p', 'f'][-magnitude-1]
        else:
            return f'$10^{magnitude*3}\\cdot$'
        
    def printTable(self):
        print(f"title: {self.title}, xlabel: {self.xlabel}, ylabel: {self.ylabel}")
        print(f"X: {self.x}, x_err: {self.x_err}.\nY: {self.y}, y_err: {self.y_err}")
    
    
    def saveAsPDF(self, path=None, filename=None, height=0, standalone=False, colnum=4):
        if filename==None:
            #sanitized_name = sanitize_filename(self.title)
            sanitized_name = self.title
            '''
            for s in r' $=':
                sanitized_name = sanitized_name.replace(s, '_')
            '''
            filename = f"{sanitized_name}_table.tex"

        if path==None:
            path = os.path.dirname(os.path.realpath(__file__))+os.path.sep
        
        if height==0: # Quadratic 
            height = int(np.ceil(self.x.shape[0]/2))
        if height == -1: # Rectangular
            height = int(self.x.shape[0])
        width = int(np.ceil(self.x.shape[0]/height))
        
        string = ''
        if standalone:
            string = r"\documentclass{article}"+"\n"+r"\usepackage{multirow}"+"\n"+r"\usepackage{geometry}"+"\n"
            string += r"\geometry{left=0cm, right=0cm, top=0cm, bottom=0cm}"+"\n"
            string += r"\begin{document}"+"\n"
        string += r"\begin{tabular}{"
        for i in range(width):
            string += r"|p{"+f"%0.2f "%(width/2/colnum-0.05)+r"\textwidth}|p{"+f"%0.2f "%(width/2/colnum-0.05)+r"\textwidth}|"
        string += r"}"
        string +="\n"+r"\hline"+"\n"+r"\multicolumn{"+f"{int(width*2)}"+r"}{|c|}{"+f"{self.title}"+r"}"+r"\\"+"\n"+r"\hline"+"\n"
        for i in range(width):
            string += f"{self.x_unit_scaled}&{self.y_unit_scaled}"
            if i<width-1:
                string += "&"
        string += r"\\"+"\n"+r"\hline"+"\n"
        x, xerr = roundToError(self.x/self.x_scaling, self.x_err/self.x_scaling)
        y, yerr = roundToError(self.y/self.y_scaling, self.y_err/self.y_scaling)
        for i in range(height):
            for k in range(width):
                if k*height+i<self.x.shape[0]:
                    string += r"$\num{"+f"{x[k*height+i]}"+r"\pm "+f"{xerr[k*height+i]}"+r"}$"+r"&$\num{"+f"{y[k*height+i]}"+r"\pm "+f"{yerr[k*height+i]}"+r"}$"
                    if(k<width-1):
                        string += "&"
                else:
                    string += "&"
            string+=r"\\"+"\n"
        string += r"\hline"+"\n"+r"\end{tabular}"
        if standalone:
            string+="\n"+r"\end{document}"
        with open(path+filename, "w") as file:
            file.write(string)
        
    def getScaling(self):
        return 10**(np.floor(np.log10(np.abs(self.x).max()))//3*3), 10**(np.floor(np.log10(np.abs(self.y).max()))//3*3)

def CSV_to_PDF(pandas_table, output_file_name):
    html_table = pandas_table.to_html()

    options = {'page-size': 'Letter', 'margin-top': '0mm',
               'margin-right': '0mm', 'margin-bottom': '0mm', 'margin-left': '0mm'}

    pdfkit.from_string(html_table, output_file_name, options=options)


def setSpace(axis, tables, border=0.1, resolution=0.1, path=None):
    if type(tables) == list:
        temp_X = np.array([])
        temp_Y = np.array([])
        temp_Xerr = np.array([])
        temp_Yerr = np.array([])
        for old_table in tables:
            temp_X = np.append(temp_X, old_table.x)
            temp_Y = np.append(temp_Y, old_table.y)
            temp_Xerr = np.append(temp_Xerr, old_table.x_err)
            temp_Yerr = np.append(temp_Yerr, old_table.y_err)
        table = Table(temp_X, temp_Xerr, temp_Y, temp_Yerr, xlabel=tables[0].xlabel, ylabel=tables[0].ylabel)
    elif type(tables) == Table:
        table = tables
    else:
        raise ValueError('Wrong type of tables: only list of olib.Table or olib.Table. Got '+f'{type(tables)}.')

        
    axis.set_xlabel(table.x_label_scaled)
    axis.set_ylabel(table.y_label_scaled)
    axis.set_title(table.title)
    
    X = table.x/table.x_scaling
    Y = table.y/table.y_scaling
    
    xRange = np.max(X)-np.min(X)
    yRange = np.max(Y)-np.min(Y)

    xstep = np.round((np.abs(xRange) * resolution)*2, -1)/2
    ystep = np.round((np.abs(yRange) * resolution)*2, -1)/2

    xLimits = np.array([math.floor((np.min(X)-xRange*border)/10)*10,
               math.ceil((np.max(X)+xRange*border+xstep)/10)*10])
    yLimits = np.array([math.floor((np.min(Y)-yRange*border)/10)*10,
               math.ceil((np.max(Y)+yRange*border+ystep)/10)*10])
    axis.set_xticks(np.arange(xLimits[0], xLimits[1], xstep))
    axis.set_xticks(np.arange(xLimits[0], xLimits[1], xstep * 0.2), minor=True)
    axis.set_yticks(np.arange(yLimits[0], yLimits[1], ystep))
    axis.set_yticks(np.arange(yLimits[0], yLimits[1], ystep * 0.2), minor=True)

    axis.set_xlim(xLimits)
    axis.set_ylim(yLimits)

    axis.grid(which="both")
    axis.grid(which="minor", alpha=0.4)
    axis.grid(which="major", alpha=0.7)
    
    if path!=None:
        if path=="":
            path = os.path.dirname(os.path.realpath(__file__))+os.path.sep+f"{sanitize_filename(table.title)}_plot.pdf"
        plt.savefig(path, dpi=500)

    return axis


def plotLine(axis, x, y, label=None, linewidth=0.7, linestyle="solid", color="r"):
    axis.plot(x, y, linewidth=linewidth, label=label,
              linestyle=linestyle, color=color)

    return axis


def plotData(axis, table, label=None, capsize=4, elinewidth=1.5, fmt="d", model=None, p0=None, color='black', fit_color='red', xscaleing=1, yscaleing=1, errorbar=True, print_data=True, data_path="data.txt", show_band=True, connect=False):
    if print_data:
        pass
    plot = None
    model_returns = None

    if model:
        model_returns = fit_function(model, table, p0=p0, plot=True, color=fit_color, show_band=show_band)
        plt.legend(facecolor="white", framealpha=1, fancybox=False)

    if errorbar:
        axis.errorbar(table.x/table.x_scaling, table.y/table.y_scaling, table.y_err/table.y_scaling, table.x_err/table.x_scaling, label=label, capsize=capsize, elinewidth=elinewidth, fmt=fmt, color=color)
    else:
        axis.scatter(table.x/table.x_scaling, table.y/table.y_scaling, label=label, marker="x", color=color, linewidth=0.6)

    if connect:
        print('<<<<<<<<<<<< connect')
        axis.plot(table.x/table.x_scaling, table.y/table.y_scaling, color=color, linewidth=0.6)

    if label != None:
        plt.legend(facecolor="white", framealpha=1, fancybox=False)

    axis.set_xlabel(table.xlabel + ' / ' + table.x_unit_scaled)
    axis.set_ylabel(table.ylabel + ' / ' + table.y_unit_scaled)
    axis.set_title(table.title)
    axis.grid()
    
    return axis if model_returns == None else [axis, model_returns]

    
    
