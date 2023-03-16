#!/usr/bin/env python3

import json
import logging
import numpy
import re
import yaml


logging.basicConfig(format = "[%(levelname)s] [%(asctime)s] %(message)s", level = logging.INFO)
logger = logging.getLogger("mylogger")


# Factorize to a utils file later
def load_config(cfgfile) :
    
    #if (os.path.isfile(cfgfile)) :
    
    with open(cfgfile, "r") as fopen :
        
        content = fopen.read()
        
        if (cfgfile.endswith(".yml") or cfgfile.endswith(".yaml")) :
            
            logger.info(f"Loading yaml config: {cfgfile}")
            d_loadcfg = yaml.load(content, Loader = yaml.FullLoader)
        
        elif (cfgfile.endswith(".json")) :
            
            logger.info(f"Loading json config: {cfgfile}")
            d_loadcfg = json.loads(content)
        
        else :
            
            logger.error(f"Invalid config provided: {cfgfile}")
            exit(1)
    
    return d_loadcfg


def parse_stau_samplestring(s) :
    
    rgx = re.compile("stau(\d+)_lsp(\d+)_ctau(\w+)")
    
    mstau, mlsp, ctau = rgx.findall(s)[0]
    
    result = {
        "mstau": float(mstau),
        "mlsp": float(mlsp),
        "ctau": ctau,
    }
    
    return result


def get_stau_xsec(samplestr, xsecfile) :
    
    arr_xsec = numpy.loadtxt(xsecfile, delimiter = ",", converters = {1: eval})
    
    d_stau_param = parse_stau_samplestring(samplestr)
    mstau = d_stau_param["mstau"]
    
    xsec = None
    
    # 1st column is the stau mass
    rowidx = numpy.where(arr_xsec[:, 0] == mstau)[0]
    
    # Empty, i.e. mstau not present
    if (not len(rowidx)) :
        
        logger.error(f"mstau {mstau} not found in {xsecfile}")
        exit(1)
    
    rowidx = rowidx[0]
    
    # 2nd colum is the xsec
    xsec = arr_xsec[rowidx, 1]
    
    return xsec


def get_stau_xsec_dict(l_samplestr, xsecfile) :
    
    d_xsec = {}
    
    for samplestr in l_samplestr :
        
        d_xsec[samplestr] = get_stau_xsec(samplestr, xsecfile)
    
    return d_xsec