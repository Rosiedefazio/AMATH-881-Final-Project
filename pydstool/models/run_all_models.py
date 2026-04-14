import os
import runpy
from pathlib import Path

import matplotlib


matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))

import matplotlib.pyplot as plt


MODEL_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = Path(os.environ.get('PYDSTOOL_MODEL_PLOT_DIR', MODEL_DIR / 'output'))
RUNNER_NAME = Path(__file__).name


def save_current_plot(output_path):
    figure_numbers = plt.get_fignums()
    if not figure_numbers:
        raise RuntimeError('No matplotlib figure was created.')
    plt.figure(figure_numbers[-1])
    plt.savefig(output_path, dpi=160, bbox_inches='tight')
    print('Saved plot to %s' % output_path)


def run_model(script_path):
    output_path = OUTPUT_DIR / ('%s.png' % script_path.stem)
    original_show = plt.show
    previous_savefig = os.environ.get('PYDSTOOL_SAVEFIG')
    did_save = {'value': False}

    def save_instead_of_show(*args, **kwargs):
        if not did_save['value']:
            save_current_plot(output_path)
            did_save['value'] = True
        plt.close('all')

    print('Running %s...' % script_path.name)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.close('all')
    os.environ['PYDSTOOL_SAVEFIG'] = str(output_path)
    plt.show = save_instead_of_show

    try:
        runpy.run_path(str(script_path), run_name='__main__')
        if not did_save['value'] and not output_path.exists():
            save_current_plot(output_path)
    finally:
        plt.show = original_show
        plt.close('all')
        if previous_savefig is None:
            os.environ.pop('PYDSTOOL_SAVEFIG', None)
        else:
            os.environ['PYDSTOOL_SAVEFIG'] = previous_savefig


def main():
    scripts = sorted(
        path for path in MODEL_DIR.glob('*.py')
        if path.name != RUNNER_NAME and not path.name.startswith('_')
    )
    if not scripts:
        raise RuntimeError('No model scripts were found under %s' % MODEL_DIR)

    for script_path in scripts:
        run_model(script_path)

    print('Generated %d plot(s) in %s' % (len(scripts), OUTPUT_DIR))


if __name__ == '__main__':
    main()
