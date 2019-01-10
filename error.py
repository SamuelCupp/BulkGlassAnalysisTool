class Error(Exception):
    """Base class for exceptions in this module."""
    pass
        
class InputError(Error):
    """
    Exception for errors raised during input phase.
    This includes reading input files and constructing
    the initial network.
    """
    def __init__(self,message):
        self.message = message

class PackageError(Error):
    """
    Exception for errors raised while interacting with
    chemistry packages.
    """
    def __init__(self,message):
        self.message = message
