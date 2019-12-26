def chr_formater(chr):
    '''
    If chromosome number cannot be directly cast to an int then the formater
    attempts to put into an int format. Currently will work if the chr num
    is the last thing in the string. If cannot format returns 0.
    '''
    try:
        return int(chr)
    except ValueError:
        for i, _ in enumerate(chr):
            try:
                return int(chr[i:])
            except ValueError:
                continue
    return 0


class feature():

    def __init__(self, ID, chr, position, type='None'):
        self.ID = ID
        self.chr = chr_formater(chr)
        self.position = int(position)
        self.type = type

    def __eq__(self, other):
        if isinstance(other, feature) and \
           self.ID == other.ID and \
           self.chr == other.chr and \
           self.position == other.position and \
           self.type == other.type:
            return True
        else:
            return False

    def __ne__(self, other):  # not needed for python3
        return not self.__eq__(other)
