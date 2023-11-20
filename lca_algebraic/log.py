import logging
import os

LOGLEVEL = os.environ.get('LOGLEVEL', 'INFO').upper()
logging.basicConfig(level=LOGLEVEL, format="[%(levelname)s] %(message)s")
logger = logging.getLogger()

