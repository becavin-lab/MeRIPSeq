#!/usr/bin/env python

import configparser

PATH='/Users/cbecavin/Documents/m6aAkker/'
config = configparser.RawConfigParser()
config.add_section('EnvVariables')
config.set('EnvVariables','PATH',PATH)
config.set('EnvVariables','cutOffMacs',10)
config.set('EnvVariables','cutOffMacs2',5)
config.set('EnvVariables','rawdataPath',PATH+'RNASEQ_raw/')
config.set('EnvVariables','fastqpath',PATH+'FastQ/')
config.set('EnvVariables','rawfastqpath',PATH+'FastQ_raw/')
config.set('EnvVariables','macsPath',PATH+'PeakDetection/Seq3_MACS2/')
config.set('EnvVariables','mappingPath',PATH+'Mapping/')
config.set('EnvVariables','expressionPath',PATH+'Expression/M6a_Input/tables/')
config.set('EnvVariables','homerPath','/usr/local/bioinf/homer/bin/')

config.add_section('ExpDesign')
config.set('ExpDesign', 'summaryFile', PATH+'Full_exp_design.csv')
config.set('ExpDesign', 'groups', ['AkkerGermFree_ZT11','AkkerGermFree_ZT3','GermFree_ZT11','GermFree_ZT3','ConvMouse_ZT11','ConvMouse_ZT3'])
config.set('ExpDesign', 'bioconds', {'AkkerGermFree','GermFree','ConvMouse'})
config.set('ExpDesign', 'comparison1', {'AkkerGermFree_ZT11':'GermFree_ZT11','AkkerGermFree_ZT11':'ConvMouse_ZT11','GermFree_ZT11':'ConvMouse_ZT11'})
config.set('ExpDesign', 'comparison2', {'AkkerGermFree_ZT3':'GermFree_ZT3','AkkerGermFree_ZT3':'ConvMouse_ZT3','GermFree_ZT3':'ConvMouse_ZT3'})

# Writing our configuration file to 'example.cfg'
with open('m6aProject.cfg', 'w') as configfile:
    config.write(configfile)
    
with open('m6aProject.cfg', 'r') as configfile:
    configNew = configparser.ConfigParser()
    configNew.readfp(configfile)
    configNew.read(configfile)
    print(configNew.sections())
    print(configNew.get('EnvVariables','PATH'))