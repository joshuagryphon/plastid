#!/usr/bin/env python

class UniqueFIFO(object):
    """FIFO of unique objects. In other words, if an element already present
    in the FIFO is appended to the FIFO, it is moved to the right end and 
    no element is removed from the FIFO. Elements are only removed when a 
    element not present in the FIFO is appended to the right end, and when
    the number of elements in the FIFO exceeds `self.max_size`.

        
    Parameters
    ----------
    idx : int
        Index of item to fetch
    
    Returns
    -------
    object
        Nth item in the |UniqueFIFO|
        
    
    Attributes
    ----------
    max_size : int
        Maximum size of |UniqueFIFO|
    """
    def __init__(self,size):
        assert size > 0
        self.max_size = size
        self._elements = []
        
    def __getitem__(self,idx):
        """Fetch Nth item in the |UniqueFIFO|
        
        Parameters
        ----------
        idx : int
            Index of item to fetch
        
        Returns
        -------
        object
            Nth item in the |UniqueFIFO|
        """
        return self._elements[idx]
    
    def __contains__(self,el):
        """Determine whether an element is in the |UniqueFIFO|
        
        Parameters
        ----------
        el : hashable object
            Query object
        
        Returns
        -------
        bool
            `True` if `el` is in the |UniqueFIFO|, `False` otherwise
        """
        return el in self._elements
    
    def __iter__(self):
        """Iterate over elements in the |UniqueFIFO|"""
        return iter(self._elements)
    
    def __len__(self):
        """Return number of objects presently in the |UniqueFIFO|
        
        Returns
        -------
        int
            Number of elements presently in the |UniqueFIFO|, must be less than
            or equal to `self.max_size`
        """
        return len(self._elements)
    
    def append(self,el):
        """Append an item to the |UniqueFIFO|. If the item is already present
        in the |UniqueFIFO|, it is moved to the right end of the FIFO, and the
        length of the |UniqueFIFO| is unchanged. Otherwise, the element is
        appended to the right end. If the length of the |UniqueFIFO| then 
        exceeds `self.max_size`, the fist element is popped.
        
        Parameters
        ----------
        el : hashable
            Any hashable object
        """
        if el in self:
            # move to front if present
            n = self._elements.index(el)
            self._elements = self._elements[:n] + self._elements[n+1:] + [el]
        else:
            # otherwise append, removing first element
            # if we are too long
            if len(self._elements) >= self.max_size:
                self._elements = self._elements[1:]
                self._elements.append(el)
            else:
                self._elements.append(el)
    
    def __str__(self):
        return str(self._elements)

    def __repr__(self):
        return "<UniqueFifo %s>" % repr(self._elements)