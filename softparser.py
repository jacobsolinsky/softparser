#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 27 11:07:33 2019

@author: Lab
"""
from itertools import chain
from warnings import warn
import re
import pandas as pd
import numpy as np
from ftplib import FTP
from io import StringIO
import tempfile
import gzip
from bokeh.plotting import figure, curdoc
from bokeh.io import show, output_notebook
from functools import lru_cache
from bokeh.models import widgets

class AttributeSet:
    def __init__(self, obligations, flags, empty_lists, full_lists):
        self.obligations, self.flags, self.empty_lists, self.full_lists = \
        obligations, flags, empty_lists, full_lists
        self.dict = {}
        self.dict['data_table_header'] = {}
        self.dict['rows'] = []
        self.dict['data_table'] = None
        self.dict['has_data_table'] = False
        for i in chain(obligations, flags):
            self.dict[i] = None
        for i in chain(empty_lists, full_lists):
            self.dict[i] = []
    
    def __setitem__(self, key, value):
        if key not in self.dict:
            self.dict[key] = []
        if type(self.dict[key]) == list:
            self.dict[key].append(value)
        elif self.dict[key] is None:
            self.dict[key] = value
        elif key in chain(self.obligations, self.flags):
            warn(f"Multiple values for {key}, there should be only one")

    def __getitem__(self, key):
        return self.dict[key]
    
    def check(self):
        for i in self.obligations:
            if self.dict[i] is None:
                warn(f"No value for {i} found")
        for i in self.full_lists:
            if self.dict[i] == []:
                warn(f"no values for {i} found")
    
    
    def __repr__(self):
        return self.dict.__repr__()


platform_attribute_set = {
       "obligations":["Platform_title",
                      "Platform_distribution",
                      "Platform_technology",
                      "Platform_manufacturer",
                      "platform_table_begin",
                      "platform_table_end"],
       "flags":["Platform_support",
                "Platform_coating",
                "Platform_geo_accession",],
       "empty_lists":["Platform_catalog_number",
                      "Platform_web_link",
                      "Platform_description",
                      "Platform_contributor",
                      "Platform_pubmed_id",
                      ],
       "full_lists":["Platform_organism", 
                     "Platform_manufacture_protocol",],
       }
null_attribute_set = {"obligations":[], "flags":[], "empty_lists":[], "full_lists":[]}


class SoftFile:
    label_value_regex = r"([^=]+) = ([^\n]*)?"
    label_novalue_regex = r"[^\n]+"
    
    def entity_indicator_line(self, line):
        groups = re.match(self.label_value_regex, line)
        if groups[1] == "PLATFORM":
            self.entity_dict[(groups[1], groups[2])] = AttributeSet(**platform_attribute_set)
        else:
            self.entity_dict[(groups[1], groups[2])] = AttributeSet(**null_attribute_set)
        self.current_entity = (groups[1], groups[2])

    def entity_attribute_line(self, line):
        groups = re.match(self.label_value_regex, line)
        if not groups:
            groups = re.match(self.label_novalue_regex, line)
            self.entity_dict[self.current_entity][groups[0]] = ""
            if groups[0] == "platform_table_begin":
                self.entity_dict[self.current_entity]["has_data_table"] = True
            else:
                self.construct_data_frame()
        else:
            self.entity_dict[self.current_entity][groups[1]] = groups[2]

    def data_table_header_description_line(self, line):
        groups = re.match(self.label_value_regex, line)
        self.entity_dict[self.current_entity]["data_table_header"][groups[1]] = groups[2]
    
    def other_line(self, line):
        self.entity_dict[self.current_entity]["rows"].append(line)

    def __init__(self, accession, full=False):
        self.header, self.has_data_table = False, False
        accession = accession.upper()
        shortpath = re.sub("\d{1,3}$", "nnn", accession)
        path = {
                "GDS": f"geo/datasets/{shortpath}/{accession}/soft/{accession}{'_full' if full else ''}.soft.gz",
                "GPL": f"geo/platforms/{shortpath}/{accession}/soft/{accession}_family.soft.gz",
                "GSE": f"geo/series/{shortpath}/{accession}/soft/{accession}_family.soft.gz",
                }[accession[:3]]
        self.entity_dict = {}
        self.ftp = FTP('ftp.ncbi.nlm.nih.gov')
        self.ftp.login()
        zipfile = tempfile.TemporaryFile()
        with zipfile:
            self.ftp.retrbinary(f'RETR {path}', zipfile.write)
            zipfile.seek(0)
            with gzip.open(zipfile, 'rt') as g:
                for line in g:
                    self.lineclassify(line)
        self.ftp.close()
        for v in self.entity_dict.values():
            v.check()

    def lineclassify(self, line):
        line_mappings = {
            "^": self.entity_indicator_line,
            "!": self.entity_attribute_line,
            "#": self.data_table_header_description_line,
            }
        if line[0] == "^":
            print(line)
        line_mappings.get(line[0], self.other_line)(line[1:])

    def construct_data_frame(self):
        try:
            a = self.entity_dict[self.current_entity]["rows"]
            if a:
                self.entity_dict[self.current_entity]["data_table"] = pd.read_table(
                        StringIO(''.join(a)))
        except pd.io.common.EmptyDataError:
            self.entity_dict[self.current_entity]["data_table"] = None
            self.has_data_table= False
 
    
    def get_entity_of_type(self, type_):
        retval = {}
        for k, v in self.entity_dict.items():
            if k[0] == type_:
                retval[k] = v
        return retval
    
    @property
    def platforms(self):
       return self.get_entity_of_type("PLATFORM")
   
    @property
    def series(self):
        return self.get_entity_of_type("SERIES")

    @property
    def samples(self):
        return self.get_entity_of_type("SAMPLE")
    
    @property
    @lru_cache(maxsize = 1)
    def rank_normalized_expression(self):
        retval = {}
        for key, sample in self.samples.items():
            retval[key[1]] = sample['data_table']['VALUE'].rank()
        return pd.DataFrame(retval)
    
    def plot_rank_normalized(self, row):
        curdoc().clear()
        s = widgets.Select(title="Gene:", options = list(next(iter(thing.samples.values()))\
                                                         ['data_table']['D_REF']))
        p = figure(plot_width = 600, plot_height = 600, 
           title = 'Example Glyphs',
           x_axis_label = 'X', y_axis_label = 'Y')
        y = self.rank_normalized_expression.loc[row, :]
        p.square(range(len(y)), y)
        show(p, 'firefox')
        show(s, 'firefox')
    
    