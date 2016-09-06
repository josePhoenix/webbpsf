import os.path

__all__ = ['PRD_VERSION', 'DATA_ROOT', 'PRD_DATA_ROOT']

DATA_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
PRD_VERSION = 'PRDDEVSOC-D-012'  # updated 2016-04-13
PRD_DATA_ROOT = os.path.join(DATA_ROOT, PRD_VERSION)
