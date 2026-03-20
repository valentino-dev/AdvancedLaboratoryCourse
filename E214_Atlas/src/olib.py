#from pathvalidate import sanitize_filename
import os
#import pdfkit
import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def roundToError(x, xerr=[""], digits=2):
    #if xerr == [""]:
        #xerr = x
    #print(xerr.shape)
    error_digits = np.zeros_like(xerr).astype("int")
    #print("esr", xerr, xerr!=0)
    #print(error_digits[xerr!=0])
    error_digits[xerr!=0] = np.floor(np.log10(xerr[xerr!=0])).astype("int")
    #error_digits[error_digits**2==0] = 0
    roundedx = []
    roundedxerr = []
    #print()
    for i in range(x.shape[0]):
        roundedx.append(np.round(x[i], digits-error_digits[i]-1))
    for i in range(xerr.shape[0]):
        roundedxerr.append(np.round(xerr[i], digits-error_digits[i]-1))
    return np.array(roundedx), np.array(roundedxerr)

class model:
    def __init__(self) -> None:
        pass

    def fit(self, X, Y):
        pass

    def predict(self, x) -> np.ndarray:
        pass


class Table:
    def __init__(self, X, Xerr, Y, Yerr, title="title", xlabel="x", ylabel="y"):
        self.X = X
        self.Xerr = Xerr
        self.Y = Y
        self.Yerr = Yerr
        
        self.title = title
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.x_scaling, self.y_scaling = self.getScaling()
        self.x_label_scaled, self.y_label_scaled = self.getScaledLables()
        
    def printTable(self):
        print(f"title: {self.title}, xlabel: {self.xlabel}, ylabel: {self.ylabel}")
        print(f"X: {self.X}, Xerr: {self.Xerr}.\nY: {self.Y}, Yerr: {self.Yerr}")
    
    
    def saveAsPDF(self, path=None, filename=None, height=0, standalone=False, colnum=2):
        if filename==None:
            sanitized_name = sanitize_filename(self.title)
            for s in r' $=':
                sanitized_name = sanitized_name.replace(s, '_')
            filename = f"{sanitized_name}_table.tex"

        if path==None:
            path = os.path.dirname(os.path.realpath(__file__))+os.path.sep
        
        if height==0: # Quadratic 
            height = int(np.ceil(self.X.shape[0]/2))
        if height == -1: # Rectangular
            height = int(self.X.shape[0])
        width = int(np.ceil(self.X.shape[0]/height))
        
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
            string += f"{self.x_label_scaled}&{self.y_label_scaled}"
            if i<width-1:
                string += "&"
        string += r"\\"+"\n"+r"\hline"+"\n"
        x, xerr = roundToError(self.X*self.x_scaling, self.Xerr*self.x_scaling)
        y, yerr = roundToError(self.Y*self.y_scaling, self.Yerr*self.y_scaling)
        for i in range(height):
            for k in range(width):
                if k*height+i<self.X.shape[0]:
                    string += r"$\num{"+f"{x[k*height+i]}"+r"\pm "+f"{xerr[k*height+i]}"+r"}{}$"+r"&$\SI{"+f"{y[k*height+i]}"+r"\pm "+f"{yerr[k*height+i]}"+r"}$"
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
        return 10**(-np.floor(np.log10(np.abs(self.X).max()))+2), 10**(-np.floor(np.log10(np.abs(self.Y).max()))+2)
    
    def getScaledLables(self):
        return self.xlabel+r"$\cdot 10^{"+f"{int(np.log10(self.x_scaling))}"+r"}$", self.ylabel+r"$\cdot 10^{"+f"{int(np.log10(self.y_scaling))}"+r"}$"
        




class LinearRegression(model):
    def __init__(self):
        self.V_m = None
        self.V_n = None
        self.m = None
        self.n = None

        return None

    def weighted_average(self, Z, Zerr):
        return (Z/Zerr**2).sum()/(1/Zerr**2).sum()

    def fit(self, X, Y, Yerr):
        self.y_wa = self.weighted_average(Y, Yerr)
        self.xy_wa = self.weighted_average(X*Y, Yerr)
        self.x_wa = self.weighted_average(X, Yerr)
        self.x_sq_wa = self.weighted_average(X**2, Yerr)

        self.sigma_wa = Yerr.size/(1/Yerr**2).sum()
        self.V_m = self.sigma_wa/(X.size*(self.x_sq_wa-self.x_wa**2))
        self.V_n = self.sigma_wa / \
            (X.size*(self.x_sq_wa-self.x_wa**2))*self.x_sq_wa
        self.m = (self.xy_wa-self.x_wa*self.y_wa)/(self.x_sq_wa-self.x_wa**2)
        self.n = (self.x_sq_wa*self.y_wa-self.x_wa*self.xy_wa) / \
            (self.x_sq_wa-self.x_wa**2)

        return self

    def setParameters(self, m, V_m, n, V_n):
        self.m = m
        self.n = n
        self.V_m = V_m
        self.V_n = V_n

    def getParameters(self):
        return {"m": self.m, "V_m": self.V_m, "n": self.n, "V_n": self.V_n}
    
    def printParameter(self):
        m_error_digits = math.floor(math.log10(self.V_m**(1/2)))
        n_error_digits = math.floor(math.log10(self.V_n**(1/2)))
        print(f"m: {np.round(self.m, 1-m_error_digits)}+-{np.round(self.V_m**(1/2), 1-m_error_digits)}; n: {np.round(self.n, 1-n_error_digits)}+-{np.round(self.V_n**(1/2), 1-n_error_digits)}")
        print(f"m: {self.m}+-{self.V_m**(1/2)}; n: {self.n}+-{self.V_n**(1/2)}")

    def predict(self, x):
        return x*self.m+self.n

    def predict_most_and_least(self, x):
        return [x*(self.m+self.V_m**(1/2))+self.n+self.V_n**(1/2), x*(self.m-self.V_m**(1/2))+self.n-self.V_n**(1/2)]


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
            temp_X = np.append(temp_X, old_table.X)
            temp_Y = np.append(temp_Y, old_table.Y)
            temp_Xerr = np.append(temp_Xerr, old_table.Xerr)
            temp_Yerr = np.append(temp_Yerr, old_table.Yerr)
        table = Table(temp_X, temp_Xerr, temp_Y, temp_Yerr, xlabel=tables[0].xlabel, ylabel=tables[0].ylabel)
    elif type(tables) == Table:
        table = tables
    else:
        raise ValueError('Wrong type of tables: only list of olib.Table or olib.Table. Got '+f'{type(tables)}.')

        
    axis.set_xlabel(table.x_label_scaled)
    axis.set_ylabel(table.y_label_scaled)
    axis.set_title(table.title)
    
    X = table.X*table.x_scaling
    Y = table.Y*table.y_scaling
    
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


def plotData(axis, table, label=None, capsize=1, elinewidth=0.5, fmt=",", polyfit=0, color=None, xscaleing=1, yscaleing=1, errorbar=True, print_data=True, data_path="data.txt"):
    if print_data:
        pass
    plot = None
    if errorbar:
        axis.errorbar(table.X*table.x_scaling, table.Y*table.y_scaling, table.Yerr*table.y_scaling, table.Xerr*table.x_scaling, label=label, capsize=capsize,
                      elinewidth=elinewidth, fmt=fmt, color=color)
    else:
        axis.scatter(table.X*table.x_scaling, table.Y*table.y_scaling, label=label, marker="x", color=color, linewidth=0.6)
    plot_color = axis.get_lines()[-1].get_color()
    #print(plot_color)

    model = None

    if polyfit != 0:
        xrange = np.max(table.X)-np.min(table.X)
        x = np.linspace(np.min(table.X)-xrange*0.1, np.max(table.X)+xrange*0.1, 100)

        if polyfit == 1:
            model = LinearRegression().fit(table.X, table.Y, table.Yerr)
            y_most_and_least = model.predict_most_and_least(x)
            plotLine(axis, x*table.x_scaling, y_most_and_least[0]*table.y_scaling,
                     linestyle="dashdot", color=plot_color, linewidth=1)
            plotLine(axis, x*table.x_scaling, y_most_and_least[1]*table.y_scaling,
                     linestyle="dashdot", color=plot_color, linewidth=1)
            y = model.predict(x)
        elif polyfit > 1:
            y = np.poly1d(np.polyfit(table.X*table.x_scaling, table.Y*table.y_scaling, polyfit))(x)

        plotLine(axis, x*table.x_scaling, y*table.y_scaling, label=None, color=plot_color, linewidth=1)
        return axis, model



    return axis

    
    
