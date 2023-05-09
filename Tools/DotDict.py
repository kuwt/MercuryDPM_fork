class DotDict(dict):

    """
    Including this class into a file will enable the user to transform dict data types into a dotdict data type
    where you can access dictionaries by dot (.) notation, e.g. instead of data['time'] we use data.time.

    myDotDict = DotDict(mydict)
    """

    # __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
    __dir__ = dict.keys

    def __getattr__(*args):
        val = dict.get(*args)
        if type(val) is dict:
            return DotDict(val)
        else:
            return val