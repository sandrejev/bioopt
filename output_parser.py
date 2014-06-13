import subprocess
import os
import pandas as p
import re


class FbaParser(object):
    def __init__(self):
        self.__p = None
        self.__stdout = None
        self.status = "Infeasible"
        self.objective = float("nan")
        self.data = []

    def run(self, cmd):
        self.__p = subprocess.Popen (
            cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE
        )

        self.__stdout = self.__read_stdout(self.__p)
        s = self.__parse_status(self.__stdout)
        self.status = s['status']
        self.objective = s['objective']
        self.data = self.__parse_fluxes()

    def __parse_fluxes(self):
        s = self.__stdout
        matches = re.findall(pattern=r"(?P<name>[^\s:]+)\s+:\s+(?P<flux>[.0-9]*)\s+(?P<cost>[.0-9]*)\s+(?P<min>[.0-9]*)\s+(?P<max>[.0-9]*)", string=s)
        res = []
        for m in matches:
            mm = {'name': m[0]}
            try:
                mm['flux'] = float(m[1])
            except ValueError, e:
                mm['flux'] = float("nan")

            try:
                mm['cost'] = float(m[2])
            except ValueError, e:
                mm['cost'] = float("nan")

            try:
                mm['min'] = float(m[3])
            except ValueError, e:
                mm['min'] = float("nan")

            try:
                mm['max'] = float(m[4])
            except ValueError, e:
                mm['max'] = float("nan")

            res.append(mm)

        return res

    def __read_stdout(self, p):
        res = ""
        for line in iter(p.stdout.readline, b''):
            res += line
        p.communicate()

        return res

    def __parse_status(self, str):
        res = {'status': None, 'objective': None}
        m = re.findall(pattern="Solution status\s+=\s+([a-zA-Z]+)", string=str)
        if m:
            res['status'] = m[0]

        m = re.findall(pattern="Solution value\s+=\s+([0-9.]+)", string=str)
        res['objective'] = float(m[0]) if res['status'] == "Optimal" and m else float('nan')

        return res

    def flux(self, name):
        return [d['flux'] for d in self.data if d['name'] == name]

    def flux_max(self, name):
        return [d['max'] for d in self.data if d['name'] == name]

    def flux_min(self, name):
        return [d['min'] for d in self.data if d['name'] == name]


def flux(data, name):
    return [d[name]['flux'] if name in d else float('nan') for d in data]


def fmin(data, name):
    return [d[name]['min'] if name in d else float('nan') for d in data]


def fmax(data, name):
    return [d[name]['max'] if name in d else float('nan') for d in data]


def r_str(model, r):
    from math import isinf
    r = model.find_reaction(r)
    lb = r.bounds.lb
    ub = r.bounds.ub

    lb = '{0}inf'.format("-" if lb < 0 else "") if isinf(lb) else '{0:.5g}'.format(lb)
    ub = '{0}inf'.format("-" if ub < 0 else "") if isinf(ub) else '{0:.5g}'.format(ub)

    return '{0} [{1}, {2}]'.format(r.name, lb, ub)


class FcfrParser(object):
    def __init__(self):
        self.__p = None
        self.data = p.DataFrame()

    def run(self, cmd, debug=False):
        output = open(os.devnull, 'wb')
        if debug:
            output = subprocess.STDOUT

        self.__p = subprocess.Popen(cmd, stdout=output, stdin=output)
        self.__p.wait()

    def read(self, path):
        self.data = p.read_table("tmp/iFF708_bioopt.fcfr", header=0, index_col=0)
        self.data = self.data.applymap(self.__extract)

        return self.data

    def __extract(self, z):
        if not isinstance(z, str):
            return {'flux': float("nan"), 'min': float("nan"), 'max': float("nan"), 'status': "None"}


        z = re.search(r"([^ ]+) [\[(]([^\])]+)[\])]\s?->\s?([^;]+)", z)
        z_vals = map(FcfrParser.__tofloat, re.split(", ", z.group(2)))
        z_dict = {'flux': z_vals[0], 'min': z_vals[0], 'max': z_vals[1], 'status': z.group(3)}

        return z_dict

    @staticmethod
    def __tofloat(v):
        try:
            return float(v)
        except ValueError, e:
            return float("nan")