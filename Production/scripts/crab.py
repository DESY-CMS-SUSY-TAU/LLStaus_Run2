#!/usr/bin/env python3
from __future__ import print_function

import argparse
import subprocess
import os
import shutil
import re
import logging
import CRABClient


def configLogger(name, loglevel=logging.INFO):
    # define a Handler which writes INFO messages or higher to the sys.stderr
    logger = logging.getLogger(name)
    logger.setLevel(loglevel)
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    console.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
    logger.addHandler(console)
    logfile = logging.FileHandler('autocrab.log')
    logfile.setLevel(loglevel)
    logfile.setFormatter(logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s'))
    logger.addHandler(logfile)


logger = logging.getLogger('autocrab')
configLogger('autocrab')
_separator = '-' * 50


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def _confirm(prompt, silent_mode=False):
    if silent_mode:
        return True
    ans = raw_input('%s [yn] ' % prompt)
    if ans.lower()[0] == 'y':
        return True
    else:
        return False


def runCrabCommand(command, *args, **kwargs):
    if kwargs.get('dryrun'):
        print(command)
        return

    from CRABAPI.RawCommand import crabCommand
    try:
        return crabCommand(command, *args, **kwargs)
    except Exception as e:
        logger.error(getattr(e, 'message', repr(e)))


def parseDatasetName(dataset):
    procname, ver, tier = dataset[1:].split('/')
    ext = ''
    isMC = tier.endswith('SIM')
    if isMC:
        ver_pieces = ver.split('_')
        keep_idx = 1
        for idx, s in enumerate(ver_pieces):
            if s.startswith('mc'):
                keep_idx = idx
                break
        rlt = re.search(r'_(v[0-9]+)(_ext[0-9]+|)(_L1v[0-9]+|)(-v[0-9]+)', ver).groups()
        ext = rlt[1].replace('_', '-') + rlt[-1]
        vername = '_'.join(ver_pieces[:keep_idx]) + '_' + rlt[0] + ext
        # hack
        if 'backup' in ver:
            ext += '_backup'
        if 'new_pmx' in ver:
            ext += '_new_pmx'
    else:
        vername = ver
        ext = '_' + ver
    return procname, vername, ext, isMC


def getDatasetSiteInfo(dataset, retry=2):
    """Return dataset storage sites for given DAS query via dasgoclient"""
    import subprocess
    import time
    import json
    query = 'site dataset=%s' % dataset
    cmd = ['dasgoclient', '-query', query, '-json']
    retry_count = 0
    while True:
        logger.info('Querying DAS:\n  %s' % ' '.join(cmd) + '' if retry_count == 0 else '\n... retry %d ...' % retry_count)
        if retry_count > 0:
            time.sleep(3)
        retry_count += 1
        if retry_count > retry:
            logger.error('Failed to retrieve site info from DAS for: %s' % dataset)
            return None, None
#             raise RuntimeError('Failed to retrieve site info from DAS for: %s' % dataset)
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        outs, errs = proc.communicate()
        if errs:
            logger.error('DAS error: %s' % errs)
            continue
        else:
            on_fnal_disk = False
            sites = []
            for row in json.loads(outs):
                for rec in row.get('site', []):
                    if rec.get('kind', '') == 'Disk' and '100' in rec.get('dataset_fraction', ''):
                        site_name = rec.get('name', '')
                        if site_name:
                            if site_name.startswith('T1_US_FNAL'):
                                on_fnal_disk = True
                            elif not site_name.startswith('T1_'):
                                sites.append(str(rec.get('name', '')))
            logger.info('Found %d sites for %s:\n  %s%s' % (len(sites), dataset, ','.join(sites), ',T1_US_FNAL' if on_fnal_disk else ''))
            return on_fnal_disk, sites


def loadConfig(work_area, task_name):
    import os
    import sys
    import copy
    from importlib import import_module
    orig_path = copy.copy(sys.path)
    cfgdir = os.path.join(work_area, 'configs')
    sys.path.insert(0, cfgdir)
    m = import_module(task_name.replace('crab_', ''))
    config = m.config
    sys.path = orig_path
    return config


def writeConfig(config, work_area):
    cfgdir = os.path.join(work_area, 'configs')
    if not os.path.exists(cfgdir):
        os.makedirs(cfgdir)
    cfgpath = os.path.join(cfgdir, config.General.requestName + '.py')
    with open(cfgpath, 'w') as f:
        f.write(str(config))
    return cfgpath


def createConfig(args, dataset):
    from CRABClient.UserUtilities import config
    config = config()

    procname, vername, ext, isMC = parseDatasetName(dataset)

    config.General.requestName = procname[:100 - len(ext)] + ext
    config.General.workArea = args.work_area
    config.General.transferOutputs = True
    config.General.transferLogs = False

    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = args.pset
    config.JobType.sendExternalFolder = args.send_external
    config.JobType.allowUndistributedCMSSW = True
    config.JobType.numCores = args.num_cores
    config.JobType.maxMemoryMB = args.max_memory
    config.JobType.maxJobRuntimeMin = args.max_runtime_min
    if args.set_input_dataset:
        config.JobType.pyCfgParams = ['inputDataset=%s' % dataset]
    if len(args.input_files) > 0:
        config.JobType.inputFiles = args.input_files

    config.Data.inputDBS = 'global'
    config.Data.inputDataset = dataset
    config.Data.splitting = args.splitting
    config.Data.unitsPerJob = args.units_per_job
    if args.max_units > 0:
        config.Data.totalUnits = args.max_units
    config.Data.publication = args.publication
    config.Data.outputDatasetTag = args.tag + '_' + vername
    config.Data.allowNonValidInputDataset = True
    config.Data.outLFNDirBase = args.outputdir

    if not isMC and args.json:
        config.Data.lumiMask = args.json

    config.Site.storageSite = args.site

    options = parseOptions(args)
    if 'siteblacklist' in options:
        config.Site.blacklist = options['siteblacklist'].split(',')

    if args.fnal:
        config.Data.ignoreLocality = True
        config.Site.whitelist = ['T3_US_FNALLPC']
        config.Site.ignoreGlobalBlacklist = True

    if args.allow_remote:
        on_fnal_disk, t2_sites = getDatasetSiteInfo(dataset)
        if on_fnal_disk and len(t2_sites) < 3:
            config.General.requestName = config.General.requestName
            config.Data.ignoreLocality = True
            config.Site.whitelist = ['T2_US_*'] + t2_sites

    # write config file
    cfgpath = writeConfig(config, args.work_area)
    return config, cfgpath


def calcLumiForRecovery(config, status_dict, work_area_rsb):
    import ast
    from CRABClient.UserUtilities import getLumiListInValidFiles
    from WMCore.DataStructs.LumiList import LumiList

    cfgdir = os.path.join(work_area_rsb, 'configs')
    if not os.path.exists(cfgdir):
        os.makedirs(cfgdir)

    # get lumis of the input dataset
    lumifile = getattr(config.Data, 'lumiMask', '')
    if lumifile:
        logger.info('Lumi mask for the original dataset: %s' % lumifile)
        if lumifile.startswith('http'):
            lumiIn = LumiList(url=lumifile)
        else:
            lumiIn = LumiList(lumifile)
    else:
        logger.info('No lumi mask for the original dataset, will use the full lumi from input dataset %s' % config.Data.inputDataset)
        lumiIn = getLumiListInValidFiles(config.Data.inputDataset, dbsurl=config.Data.inputDBS)

    # get lumis of the processed dataset
    outputDataset = ast.literal_eval(status_dict['outdatasets'])[0]
    logger.info('Getting lumis in the output dataset %s' % outputDataset)
    lumiDone = getLumiListInValidFiles(outputDataset, dbsurl='phys03')
    lumiDone.writeJSON(os.path.join(cfgdir, config.General.requestName + '_lumi_processed.json'))

    outpath = os.path.abspath(os.path.join(cfgdir, config.General.requestName + '_lumiMask.json'))
    newLumiMask = lumiIn - lumiDone
    newLumiMask.writeJSON(outpath)
    return outpath


def parseOptions(args):

    def convertValue(v):
        if v.lower() == 'true':
            v = True
        elif v.lower() == 'false':
            v = False
        return v

    options = {}
    if args.options:
        prev = None
        for opt in args.options.split():
            if '=' in opt:
                k, v = opt.split('=')
                if k.startswith('--'):
                    k = k[2:]
                options[k] = convertValue(v)
            else:
                if opt.startswith('--'):
                    if prev is None:
                        prev = opt[2:]
                        continue
                    else:
                        options[prev] = True
                        prev = opt[2:]
                else:
                    options[prev] = convertValue(opt)
                    prev = None
    return options


def killjobs(args):
    import os
    for work_area in args.work_area:
        for dirname in os.listdir(work_area):
            logger.info('Kill job %s/%s' % (work_area, dirname))
            runCrabCommand('kill', dir='%s/%s' % (work_area, dirname))


def resubmit(args):
    import os
    kwargs = parseOptions(args)
    for work_area in args.work_area:
        for dirname in os.listdir(work_area):
            logger.info('Resubmitting job %s/%s with options %s' % (work_area, dirname, str(kwargs)))
            runCrabCommand('resubmit', dir='%s/%s' % (work_area, dirname), **kwargs)


def _analyze_crab_status(ret):
    # https://github.com/dmwm/CRABClient/blob/master/src/python/CRABClient/Commands/status.py
    states = {}
    statesPJ = {}  # probe jobs
    statesSJ = {}  # tail jobs
    for jobid in ret['jobs']:
        jobStatus = ret['jobs'][jobid]['State']
        if jobid.startswith('0-'):
            statesPJ[jobStatus] = statesPJ.setdefault(jobStatus, 0) + 1
        elif '-' in jobid:
            statesSJ[jobStatus] = statesSJ.setdefault(jobStatus, 0) + 1
        else:
            states[jobStatus] = states.setdefault(jobStatus, 0) + 1

    if sum(statesPJ.values()) > 0:
        if 'failed' in states and sum(statesSJ.values()) > 0:
            # do not consider failed jobs that are re-scheduled
            states.pop('failed')
        for jobStatus in statesSJ:
            states[jobStatus] = states.setdefault(jobStatus, 0) + statesSJ[jobStatus]

    return states


def status(args):
    import os
    import json
    kwargs = parseOptions(args)
    for work_area in args.work_area:
        # load status from last query
        _task_status_file = os.path.join(work_area, 'task_status.json')
        crab_task_status = {}
        if os.path.exists(_task_status_file):
            with open(_task_status_file) as f:
                crab_task_status = json.load(f)

        if args.prepare_recovery_task or args.submit_recovery_task:
            work_area_rsb = work_area.rstrip('/') + args.recovery_task_suffix
            _recovery_task_file = os.path.join(work_area_rsb, 'recovery_tasks.json')
            if not os.path.exists(work_area_rsb):
                os.makedirs(work_area_rsb)
            if args.prepare_recovery_task:
                recovery_tasks = {}
            else:
                with open(_recovery_task_file) as f:
                    recovery_tasks = json.load(f)

        jobnames = [d for d in os.listdir(work_area) if d.startswith('crab_')]
        finished = 0
        job_status = {}
        submit_failed = []
        for dirname in jobnames:
            if args.submit_recovery_task:
                if dirname not in recovery_tasks or not recovery_tasks[dirname]['resubmit']:
                    continue

            # skip jobs that are already completed
            if dirname in crab_task_status:
                if crab_task_status[dirname].get('status', '') == 'COMPLETED':
                    logger.info('Skip completed job %s' % dirname)
                    finished += 1
                    continue
            # check task status
            logger.info('Checking status of job %s' % dirname)
            ret = runCrabCommand('status', dir='%s/%s' % (work_area, dirname))
            try:
                states = _analyze_crab_status(ret)
            except:
                logger.warning('Cannot get status for job %s' % dirname)
                job_status[dirname] = '\033[1;101mUNKNOWN\033[0m'
                continue
            try:
                percent_finished = 100.*states['finished'] / sum(states.values())
            except KeyError:
                percent_finished = 0
            pcts_str = ' (\033[1;%dm%.1f%%\033[0m)' % (32 if percent_finished > 90 else 34 if percent_finished > 70 else 35 if percent_finished > 50 else 31, percent_finished)
            job_status[dirname] = ret['status'] + pcts_str + '\n    ' + str(states)
            if ret['publicationEnabled']:
                pcts_published = 100.* ret['publication'].get('done', 0) / max(sum(states.values()), 1)
                pub_pcts_str = '\033[1;%dm%.1f%%\033[0m' % (32 if pcts_published > 90 else 34 if pcts_published > 70 else 35 if pcts_published > 50 else 31, pcts_published)
                job_status[dirname] = job_status[dirname] + '\n    publication: ' + pub_pcts_str + ' ' + str(ret['publication'])
                if ret['status'] == 'COMPLETED' and pcts_published != 100:
                    ret['status'] == 'PUBLISHING'

            # save task status
            crab_task_status[dirname] = ret

            if args.prepare_recovery_task or args.submit_recovery_task:
                if args.prepare_recovery_task:
                    if ret['status'] == 'COMPLETED':
                        continue
                    # first kill the task
                    if _confirm('Kill job %s/%s and prepare a recovery task?' % (work_area, dirname), silent_mode=args.yes):
                        runCrabCommand('kill', dir='%s/%s' % (work_area, dirname))  # FIXME
                        recovery_tasks[dirname] = {'completed': percent_finished, 'resubmit': True}

                elif args.submit_recovery_task:
                    if 'KILLED' not in ret['status']:
                        skip = _confirm('Task %s/%s is not in status KILLED, wait and submit the recovery task later?' % (work_area, dirname), silent_mode=args.yes)
                        if skip:
                            continue
                    config = loadConfig(work_area, dirname)
                    config.General.workArea = work_area_rsb
                    config.Data.lumiMask = calcLumiForRecovery(config, ret, work_area_rsb)
                    cfgpath = writeConfig(config, work_area_rsb)
                    if args.dryrun:
                        print('-' * 50)
                        print(config)
                        continue
                    logger.info('Submitting recovery task for %s/%s' % (work_area, dirname))
                    cmd = 'crab submit -c {cfgpath}'.format(cfgpath=cfgpath)
                    p = subprocess.Popen(cmd, shell=True)
                    p.communicate()
                # move on to next task
                continue

            if ret['status'] == 'COMPLETED':
                finished += 1
            elif ret['dbStatus'] == 'SUBMITFAILED':
                if not args.no_resubmit:
                    logger.info('Resubmitting submit-failed job %s.' % dirname)
                    shutil.rmtree('%s/%s' % (work_area, dirname))
                    cfgpath = os.path.join(work_area, 'configs', dirname.lstrip('crab_') + '.py')
                    cmd = 'crab submit -c {cfgpath}'.format(cfgpath=cfgpath)
                    p = subprocess.Popen(cmd, shell=True)
                    p.communicate()
                    if p.returncode != 0:
                        submit_failed.append(ret['inputDataset'])
            #elif states.get('failed', 0) > 0 and 'killed' not in ret['status'].lower() and not args.no_resubmit:
            #    logger.info('Resubmitting job %s with options %s' % (dirname, str(kwargs)))
            #    runCrabCommand('resubmit', dir='%s/%s' % (work_area, dirname), **kwargs)

            if ret['publication'].get('failed', 0) > 0:
                logger.info('Resubmitting job %s for failed publication' % dirname)
                runCrabCommand('resubmit', dir='%s/%s' % (work_area, dirname), publication=True)


        logger.info('====== Summary (%s) ======\n' % (work_area) +
                     '\n'.join(['%s: %s' % (k, job_status[k]) for k in natural_sort(job_status.keys())]))
        logger.info('%d/%d jobs completed!' % (finished, len(jobnames)))
        if len(submit_failed):
            logger.warning('Submit failed:\n%s' % '\n'.join(submit_failed))

        # write job status file
        with open(_task_status_file, 'w') as f:
            json.dump(crab_task_status, f, indent=2)

        if args.prepare_recovery_task:
            with open(_recovery_task_file, 'w') as f:
                json.dump(recovery_tasks, f, indent=2)


def summary_from_log_file():
    import ast
    summary = {}
    with open('autocrab.log') as f:
        for l in f:
            if _separator in l:
                # skip all previous runs
                summary = {}
                continue
            l = l.strip()
            if l[:2] == "{'":
                s = ast.literal_eval(l)
                for k in s:
                    if k not in summary:
                        summary[k] = s[k]
                    else:
                        summary[k] += s[k]
    print(' === Summary: ', summary)


def main():

    parser = argparse.ArgumentParser('Submit crab jobs')
    parser.add_argument('-i', '--inputfile',
                        help='File with list of input datasets'
                        )
    parser.add_argument('-o', '--outputdir',
                        help='Output directory'
                        )
    parser.add_argument('-p', '--pset',
                        help='Path to the CMSSW configuration file'
                        )
    parser.add_argument('-s', '--splitting',
                        default='Automatic', choices=['Automatic', 'FileBased', 'LumiBased', 'EventAwareLumiBased'],
                        help='Job splitting method. Default: %(default)s'
                        )
    parser.add_argument('-n', '--units-per-job',
                        default=300, type=int,
                        help='Units per job. The meaning depends on the splitting. Recommended default numbers: (Automatic: 300 min, LumiBased:100, EventAwareLumiBased:100000) Default: %(default)d'
                        )
    parser.add_argument('--max-units',
                        default=-1, type=int,
                        help='Max units per job. The meaning depends on the splitting. Default: %(default)d'
                        )
    parser.add_argument('-t', '--tag',
                        default='NanoHRT',
                        help='Output dataset tag. Default: %(default)s'
                        )
    parser.add_argument('-j', '--json',
                        default='https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt',
                        help='JSON file for lumi mask. Default: %(default)s'
                        )
    parser.add_argument('--site',
                        default='T3_US_FNALLPC',
                        help='Storage site. Default: %(default)s'
                        )
    parser.add_argument('--send-external',
                        action='store_true', default=False,
                        help='Send external folder. Default: %(default)s'
                        )
    parser.add_argument('--publication',
                        action='store_true', default=False,
                        help='Publish the output dataset. Default: %(default)s'
                        )
    parser.add_argument('--set-input-dataset',
                        action='store_true', default=False,
                        help='Set the inputDataset parameter to the python config file. Default: %(default)s'
                        )
    parser.add_argument('--input-files',
                        default=[], nargs='*',
                        help='Set input files to be shipped with the CRAB jobs. Default: %(default)s'
                        )
    parser.add_argument('--work-area', nargs='*',
                        default='crab_projects',
                        help='Crab project area. Default: %(default)s'
                        )
    parser.add_argument('--num-cores',
                        default=1, type=int,
                        help='Number of CPU cores. Default: %(default)d'
                        )
    parser.add_argument('--max-memory',
                        default=2000, type=int,
                        help='Number of memory. Default: %(default)d MB'
                        )
    parser.add_argument('--max-runtime-min',
                        default=2750, type=int,
                        help='Maximum job runtime (in minutes). Values above 2750 will be automatically changed to 2750 by CRAB. Default: %(default)d min'
                        )
    parser.add_argument('--dryrun',
                        action='store_true', default=False,
                        help='Only print the commands but do not submit. Default: %(default)s'
                        )
    parser.add_argument('--fnal',
                        action='store_true', default=False,
                        help='Run at FNAL LPC. Default: %(default)s'
                        )
    parser.add_argument('--allow-remote',
                        action='store_true', default=False,
                        help='Allow jobs to run remotely under certain conditions. Default: %(default)s'
                        )
    parser.add_argument('--status',
                        action='store_true', default=False,
                        help='Check job status. Will resubmit if there are failed jobs. Default: %(default)s'
                        )
    parser.add_argument('--no-resubmit',
                        action='store_true', default=False,
                        help='Disable auto resubmit of SUBMITFAILED tasks when checking job status. Default: %(default)s'
                        )
    parser.add_argument('--resubmit',
                        action='store_true', default=False,
                        help='Resubmit jobs. Default: %(default)s'
                        )
    parser.add_argument('--kill',
                        action='store_true', default=False,
                        help='Kill jobs. Default: %(default)s'
                        )
    parser.add_argument('--prepare-recovery-task',
                        action='store_true', default=False,
                        help='Prepare recovery tasks. This will kill the current jobs. Default: %(default)s'
                        )
    parser.add_argument('--submit-recovery-task',
                        action='store_true', default=False,
                        help='Submit recovery tasks. This will check if the original job has been killed. Default: %(default)s'
                        )
    parser.add_argument('--recovery-task-suffix',
                        default='_rsb',
                        help='Suffix for the work area of the recovery tasks. Default: %(default)s'
                        )
    parser.add_argument('-y', '--yes',
                        action='store_true', default=False,
                        help='Do not ask for confirmation. Default: %(default)s'
                        )
    parser.add_argument('--options',
                        default='',
                        help='CRAB command options, space separated string. Default: %(default)s'
                        )
    parser.add_argument('--summary',
                        action='store_true', default=False,
                        help='Print job status summary from the log file. Default: %(default)s'
                        )
    args = parser.parse_args()

    if args.summary:
        summary_from_log_file()
        return

    # write a separator to distinguish between different runs
    logger.info(_separator)

    if args.status or args.prepare_recovery_task or args.submit_recovery_task:
        status(args)
        return

    if args.resubmit:
        resubmit(args)
        return

    if args.kill:
        killjobs(args)
        return

    assert(len(args.work_area) == 1)
    args.work_area = args.work_area[0]

    submit_failed = []
    request_names = {}
    with open(args.inputfile) as inputfile:
        for l in inputfile:
            l = l.strip()
            if not l or l.startswith('#'):
                continue
            dataset = [s for s in l.split() if '/MINIAOD' in s][0]
            cfg, cfgpath = createConfig(args, dataset)
            if cfg.General.requestName in request_names:
                request_names[cfg.General.requestName].append(dataset)
            else:
                request_names[cfg.General.requestName] = [dataset]
            if args.dryrun:
                print('-' * 50)
                print(cfg)
                continue
            logger.info('Submitting dataset %s' % dataset)
            cmd = 'crab submit -c {cfgpath}'.format(cfgpath=cfgpath)
            p = subprocess.Popen(cmd, shell=True)
            p.communicate()
            if p.returncode != 0:
                submit_failed.append(cfgpath)
#             runCrabCommand('submit', config=cfg)

    if len(submit_failed):
        logger.warning('Submit failed:\n%s' % '\n'.join(submit_failed))
    duplicate_names = {name:request_names[name] for name in request_names if len(request_names[name])>1}
    if len(duplicate_names):
        logger.warning('Dataset with the same request names:\n%s' % str(duplicate_names))


if __name__ == '__main__':
    main()
