import ast
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]


def shared_config(settings):
    with tempfile.TemporaryDirectory() as tmp:
        config = dict(settings, output_dir=tmp, tmp_root=tmp)
        scope = {'config': config, 'workflow': SimpleNamespace(basedir=str(ROOT))}
        source = ROOT / 'workflows/shared_config.smk'
        exec(compile(source.read_text(), str(source), 'exec'), scope)
        return scope


class HostExclusionTests(unittest.TestCase):
    def test_default_empty_and_null_preserve_original_requirements(self):
        for settings in ({}, {'lsf_exclude_hosts': []}, {'lsf_exclude_hosts': None}):
            self.assertEqual(shared_config(settings)['LSF_HOST_SELECT'], '')

    def test_single_host(self):
        self.assertEqual(shared_config({'lsf_exclude_hosts': ['r640c15']})['LSF_HOST_SELECT'],
                         "select[type==local&&hname!='r640c15'] ")

    def test_multiple_hosts_are_quoted_and_deduplicated(self):
        scope = shared_config({'lsf_exclude_hosts': ['node-1.example.org', 'r640c15', 'r640c15']})
        self.assertEqual(scope['LSF_HOST_SELECT'],
                         "select[type==local&&hname!='node-1.example.org'&&hname!='r640c15'] ")

    def test_invalid_types_and_shell_metacharacters_are_rejected(self):
        for value in ('r640c15', '', False, 0, {}, [None], [12], [''], ['a b'],
                      ['a\nb'], ["a'"], ['$(touch bad)'], ['a;exit'], ['a*']):
            with self.subTest(value=value), self.assertRaises(ValueError):
                shared_config({'lsf_exclude_hosts': value})

    def test_active_yaml_and_sample(self):
        active = yaml.safe_load((ROOT / 'yaml/config_state4.yaml').read_text())
        sample = yaml.safe_load((ROOT / 'config_sample.yaml').read_text())
        self.assertIn('r640c15', active['lsf_exclude_hosts'])
        self.assertEqual(sample['lsf_exclude_hosts'], [])

    def test_every_cluster_rule_has_host_selection(self):
        for path in sorted((ROOT / 'workflows').glob('step*.smk')):
            source = path.read_text()
            blocks = re.findall(r'^    resources:\n((?:.*\n)*?)(?=^    \w|\Z)', source, re.MULTILINE)
            for block in blocks:
                if 'queue=' in block:
                    self.assertIn('lsf_host_select=LSF_HOST_SELECT', block, str(path))

    def test_all_nested_bsub_calls_include_selection_as_one_argument(self):
        prefix = shared_config({'lsf_exclude_hosts': ['r640c15', 'node-2']})['LSF_HOST_SELECT']
        count = 0
        for path in sorted((ROOT / 'workflows').glob('step*.smk')):
            for line in path.read_text().splitlines():
                if '-R ' not in line:
                    continue
                count += 1
                self.assertIn('{resources.lsf_host_select}', line, str(path))
                rendered = line.replace('{resources.lsf_host_select}', prefix).strip().rstrip('\\').strip()
                words = shlex.split(rendered)
                self.assertEqual(len(words), 2)
                self.assertEqual(words[0], '-R')
                self.assertTrue(words[1].startswith(prefix))
                self.assertIn('span[hosts=1]', words[1])
        self.assertEqual(count, 3)

    def test_profile_and_submit_script_preserve_full_resource_expression(self):
        profile = yaml.safe_load((ROOT / 'profile/lsf/config.yaml').read_text())
        for hosts in ([], ['r640c15', 'node-2']):
            with self.subTest(hosts=hosts), tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                stub = tmp / 'bsub'
                shutil.copyfile(ROOT / 'tests/fixtures/bsub_stub.py', stub)
                stub.chmod(0o755)
                jobscript = tmp / 'job script.sh'
                jobscript.write_text('#!/bin/bash\necho test\n')
                record = tmp / 'record.json'
                prefix = shared_config({'lsf_exclude_hosts': hosts})['LSF_HOST_SELECT']
                resources = SimpleNamespace(queue='long short', mem_mb=500, cpus=1,
                                            walltime='8:00', lsf_host_select=prefix)
                command = profile['cluster-generic-submit-cmd'].format(resources=resources, rule='test', jobid=1)
                argv = shlex.split(command) + [str(jobscript)]
                env = dict(os.environ, PATH=str(tmp) + os.pathsep + os.environ['PATH'], BSUB_TEST_RECORD=str(record))
                result = subprocess.run(argv, cwd=ROOT, env=env, capture_output=True, text=True, check=True)
                self.assertEqual(result.stdout.strip(), '123456')
                saved = json.loads(record.read_text())
                expected = prefix + 'rusage[mem=500]'
                self.assertEqual(saved['argv'][saved['argv'].index('-R') + 1], expected)
                self.assertEqual(saved['argv'][saved['argv'].index('-q') + 1], 'long short')
                self.assertEqual(saved['stdin'], jobscript.read_text())
                self.assertEqual(saved['argv'][:2], ['-u', ''])
                subprocess.run(['bash', str(ROOT / 'profile/lsf/submit.sh'), '-q', 'long short',
                                '-R', expected, str(jobscript)], env=env, check=True, capture_output=True)
                saved = json.loads(record.read_text())
                self.assertEqual(saved['argv'][saved['argv'].index('-q') + 1], 'long short')

    def test_all_workflows_parse_with_installed_snakemake(self):
        from snakemake import parser
        from snakemake.sourcecache import LocalSourceFile
        workflow = SimpleNamespace(sourcecache=SimpleNamespace(open=lambda src: open(src.get_path_or_uri(secret_free=True))))
        for path in sorted((ROOT / 'workflows').glob('step*.smk')):
            code, _ = parser.parse(LocalSourceFile(str(path)), workflow, {})
            ast.parse(code)


if __name__ == '__main__':
    unittest.main()
