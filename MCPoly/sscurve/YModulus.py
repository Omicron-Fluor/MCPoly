import re
import numpy as np
import pandas as pd
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression

def YModulus(file):
    '''
    Try to calculate the Young's Modulus based on LinearRegression.
    YModulus(file)
    file: File Name
    Unit: nN
    Use the calculated .txt file from sscurve.single or sscurve.multiple to get Young's Modulus.
    '''
    datum = pd.read_csv(file)
    b = re.search(r'[A-Z]+[a-z]*[A-Z]*[a-z]*[0-9]*\-?[A-Z]*[0-9]*_', file)
    #print(b.group(0)[:-1])
    datum = datum.rename(columns = {
        'Strain Length(%)':'Strain Length ({0})'.format(b.group(0)[:-1])})
    data = datum['Strain Length ({0})'.format(
        b.group(0)[:-1])].to_numpy().reshape(-1, 1)
    target = datum['Stress Force(nN)'].to_numpy().reshape(-1, 1)
    j = 3
    while 1:
        data2 = data[1:j]
        target2 = target[1:j]
        model = LinearRegression()
        model.fit(data2, target2)
        y = model.coef_[0][0]*data2 + model.intercept_
        test_score = r2_score(target2, y)
        #print(model.coef_[0][0], model.intercept_[0], test_score)
        if test_score < 0.99:
            break
        k = model.coef_[0][0]
        b = model.intercept_[0]
        j = j + 1
    return k