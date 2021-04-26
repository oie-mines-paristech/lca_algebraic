# This file provide additionnal database handler for remote single file DB
import os

from bw2data.backends import SingleFileDatabase
from bw2data import config


class RemoteSingleFileDB(SingleFileDatabase) :
    """ Targets a single database file (pickle) in a remote folder"""

    backend = "remotesinglefile"

    def __init__(self, name, filename=None, folder=None, version=None):
        """ User should provide either full filename of the pickle file, or containing folder : in the later case,
        the file name is derived from the database name as done by SingleFileDatabase, by adding hash of the name."""

        #assert not (filename is None and folder is None), "You should provide either 'filename' : full path to pickle file " \
        #                                            "or 'folder : Path of 'intermediate' folder of remote project."

        super(RemoteSingleFileDB, self).__init__(name)
        self.register(filename=filename, folder=folder, version=version)

    def filepath_intermediate(self, version=None):
        if self.metadata.get("filename") is not None :
            return self.metadata["filename"]
        else :
            return os.path.join(
                self.metadata["folder"],
                self.filename_for_version(version) + ".pickle")

# Register it
config.backends[RemoteSingleFileDB.backend] = RemoteSingleFileDB