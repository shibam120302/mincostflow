import ctypes.util
import ctypes
import binascii
import numpy

libmcf_found = ctypes.util.find_library('mincostflow')

if not libmcf_found:
    raise Exception('libmincostflow not found in LD_LIBRARY_PATH')

libmcf = ctypes.cdll.LoadLibrary(libmcf_found)

class arcID_type(ctypes.Structure):
    '''
    Arc id is defined as a c-like struct with two variables:
    
    short_channel_id (8 bytes)
    = block_height (3 bytes) || transaction_index (3 bytes) || output_index (2 bytes)
    
    part (1 byte)
    = piece_id (7 bits) || direction (1 bit) 
    
    The piece_id is an ordinal number to identify each of the multiple pieces in which a channel is
    split in the piecewise linear approximation of the cost function.
    
    The direction bit purpose is to distinguish the channel direction.
    For a channel that connects nodes A and B, we have direction = 0 if A<B
    and direction = 1 if A>B.
    '''
    _fields_ = [('short_channel_id',ctypes.c_longlong),
                ('part',ctypes.c_ubyte)]
    def from_args(arc_id: str, arc_part: int):
        '''
        A function to construct an arcID_type from a
        short channel id and 'part' index.
        '''
        a,b,c=arc_id.split('x')
        id_num = (int(a)<<40)+(int(b)<<16)+int(c)
        return arcID_type(id_num,arc_part)
    def serialize(self):
        '''
        Returns the short channel id (as a string) and the index of the channel part of the
        piecewise linear-cost channels decomposition
        '''
        scid = self.short_channel_id 
        str_scid = str( (scid>>40) & 0xffffff ) + 'x'+ str( (scid>>16) & 0xffffff )  + 'x' + str( scid & 0xffff )
        return str_scid,self.part

class nodeID_type(ctypes.Structure):
    '''
    The node id here is represented as a 33 byte array 
    representing its public key.
    '''
    _fields_ = [('k',ctypes.c_ubyte * 33)]
    def from_str(nodeid : str):
        '''
        It returns a nodeID_type from the string that represents the node public key.
        '''
        binary = binascii.unhexlify(nodeid)
        return nodeID_type.from_buffer_copy(binary)
    def serialize(self):
        '''
        It returns the string that encodes the node public key.
        '''
        s=[]
        for i in range(33):
            s.append(self.k[i])
        return bytes(s).hex()

