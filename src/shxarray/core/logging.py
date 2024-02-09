# This file is part of the shxarray software which is licensed
# under the Apache License version 2.0 (see the LICENSE file in the main repository)
# Copyright Roelof Rietbroek (r.rietbroek@utwente.nl), 2023
#


import logging
# shxarray  wide logger
logger=logging.getLogger("shxarray")

ch = logging.StreamHandler()

# create formatter
formatter = logging.Formatter('%(name)s-%(levelname)s: %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)


def debugging():
    return logger.getEffectiveLevel() == logging.DEBUG

def setInfoLevel():
    """Set logging level for both python and c++ to INFO severity"""
    logger.setLevel(logging.INFO)

def setDebugLevel():
    """Set logging level for both python and c++ to DEBUG severity"""
    logger.setLevel(logging.DEBUG)


def setWarningLevel():
    """Set logging level for both python and c++ to WARNING severity"""
    logger.setLevel(logging.WARNING)

def setErrorLevel():
    """Set logging level for both python and c++ to WARNING severity"""
    logger.setLevel(logging.ERROR)

setInfoLevel()
