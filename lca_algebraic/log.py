import logging
import os

LOGLEVEL = os.environ.get("LOGLEVEL", "INFO").upper()
# logging.basicConfig(level=LOGLEVEL, format="[%(levelname)s] %(message)s")

logger = logging.getLogger("lca_algebraic")
logger.setLevel(LOGLEVEL)

ch = logging.StreamHandler()
ch.setLevel(LOGLEVEL)
ch.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))

logger.addHandler(ch)


def debug(*args):
    msg = " ".join(str(item) for item in args)
    logger.debug(msg)


def warn(*args):
    msg = " ".join(str(item) for item in args)
    logger.warning(msg)
