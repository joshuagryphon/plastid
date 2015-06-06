#!/usr/bin/env python
"""Tests capabilities of |ArrayTable|"""
import unittest
import os
import tempfile
import copy
import numpy
from nose.plugins.attrib import attr
from yeti.util.services.exceptions import MalformedFileError
from yeti.util.array_table import ArrayTable, equal_enough, _NUMERIC_DTYPES
from yeti.util.io.filters import NameDateWriter

printer = NameDateWriter("test")

@attr(test="unit")
class TestArrayTable(unittest.TestCase):
    """TestCase for :py:class:`yeti.util.array_table`"""
    
    @classmethod
    def setUpClass(cls):
        chars = "ABCDEFGHJKLMNOPQRSTUVWXYZabcdefghijklmopqrstuvwxyz"
        my_size = 1000
        
        # test read/write over large range of float values
        exps = numpy.arange(-50,50)
        dtmp = { str(X) : (10**X)*numpy.random.random(size=my_size) for X in exps }
        
        # add integer columns
        dtmp["int1"] = numpy.random.randint(0,high=100000,size=my_size)
        dtmp["int2"] = numpy.random.randint(0,high=100000,size=my_size)
        
        # add string columns
        dtmp["string1"] = numpy.tile(None,my_size)
        dtmp["string2"] = numpy.tile(None,my_size)
        for i in range(my_size):
            dtmp["string1"][i] = chars[numpy.random.randint(len(chars))] + chars[numpy.random.randint(len(chars))] + chars[numpy.random.randint(len(chars))] 
            dtmp["string2"][i] = chars[numpy.random.randint(len(chars))] + chars[numpy.random.randint(len(chars))] + chars[numpy.random.randint(len(chars))] + chars[numpy.random.randint(len(chars))] 
             
        # add bool columns
        dtmp["bool1"] = numpy.tile(True,my_size)
        dtmp["bool2"] = numpy.tile(False,my_size)
        dtmp["bool3"] = numpy.tile(False,my_size)
        dtmp["bool4"] = numpy.tile(True,my_size)    
        dtmp["bool1"][numpy.random.randint(0,my_size,size=5000)] = False
        dtmp["bool2"][numpy.random.randint(0,my_size,size=5000)] = True
        
        # add nan column
        dtmp["nan"] = numpy.tile(numpy.nan,my_size)
        
        cls.val_dict = dtmp
        cls.table1 = ArrayTable(dtmp)
        
    def test_creation_assumptions(self):
        """Assure we can't create an ArrayTable with uneven columns"""
        dtmp = dict(a=range(10),b=range(50),c=range(10,20))
        self.assertRaises(ValueError,ArrayTable,dtmp)
        
        dtmp = dict(a=range(10),b=[])
        self.assertRaises(ValueError,ArrayTable,dtmp)
    
    def test_keys(self):
        a = ArrayTable(dict(a=range(50),b=range(50),c=range(50)))
        self.assertSetEqual(set(a.keys()), set(list("abc")))

        a["d"] = range(50,100)
        self.assertSetEqual(set(a.keys()), set(list("abcd")))
        self.assertNotEquals(set(a.keys()), set(list("abc")))

    def test_values(self):
        self.assertTrue(set([hash(tuple(X)) for X in self.table1.values()]),
                        set([hash(tuple(X)) for X in self.val_dict.values()]))

    def test_items(self):
        d = dict(a=range(50),b=range(150,200),c=range(250,300))
        a = ArrayTable(d)
        for k,v in a.items():
            self.assertTrue((a[k] == v).all())
            self.assertTrue((a[k] == d[k]).all())

    def test_shape(self):
        a = ArrayTable(dict(a=range(50),b=range(50),c=range(50)))
        self.assertTupleEqual(a.shape(), (50,3))
        
        a["d"] = range(50,100)
        self.assertTupleEqual(a.shape(), (50,4))
    
    def test_len(self):
        for _ in range(10):
            my_len = numpy.random.randint(1,300)
            dtmp = { "a" : numpy.random.random(size=my_len),
                     "b" : numpy.random.random(size=my_len),
                     "c" : numpy.random.random(size=my_len)
                     }
            a = ArrayTable(dtmp)
            self.assertEquals(len(a),my_len)
    
    def test_getitem(self):
        for k,v in self.val_dict.items():
            self.assertTrue(equal_enough(v,self.table1[k]))

    def test_setitem(self):
        a = ArrayTable(dict(a=range(50),b=range(50),c=range(50)))
        self.assertTupleEqual(a.shape(), (50,3))
        d = range(50,100)
        old_a = a["a"]
        old_b = a["b"]
        old_c = a["c"]
        
        # add new column
        a["d"] = d
        self.assertTupleEqual(a.shape(), (50,4))
        self.assertTrue(equal_enough(a["d"],numpy.array(d)))
        
        # make sure other columns didn't change
        self.assertTrue(equal_enough(a["a"],old_a))
        self.assertTrue(equal_enough(a["b"],old_b))
        self.assertTrue(equal_enough(a["c"],old_c))
        
        # make sure we can replace existing columns
        new_a = range(40,90)
        a["a"] = new_a
        self.assertTrue(equal_enough(a["a"],numpy.array(new_a)))
        self.assertFalse(equal_enough(a["a"],numpy.array(old_a)))
     
    def test_setitem_assumptions(self):
        """Assure we can't add columns that are mismatched in length"""
        a = ArrayTable(dict(a=range(50),b=range(50),c=range(50)))
        self.assertRaises(ValueError,a.__setitem__,"d",range(51))

    def test_get_rows(self):
        keyorder = sorted(self.table1.keys())
        
        # fetch whole table --------------------------------
        rows = list(self.table1.get_rows(keyorder=keyorder))
        self.assertTrue(len(rows),len(self.table1))
        for i,row in enumerate(rows):
            # make sure rows are correct length
            self.assertEquals(len(row),len(self.table1.keys()))
            
            # make sure cells are correct
            for j, col in enumerate(keyorder):
                if not isinstance(row[j],float) or (isinstance(row[j],float) and not numpy.isnan(row[j])):
                    val1 = row[j]
                    val2 = self.table1[col][i]
                    val3 = self.val_dict[col][i]
                    self.assertEquals(val1,val2,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val2))
                    self.assertEquals(val1,val3,
                                      "Row cell value does not match Dictionary cell value (col %s): %s vs %s" % (col,val1,val3))

        # fetch masking rows -------------------------------
        mymask  = self.table1["bool1"] == True
        mymask &= self.table1["bool2"] == True
        mymask |= self.table1["int1"] > 70000
        self.assertGreater(mymask.sum(),0)
        
        rows = list(self.table1.get_rows(keyorder=keyorder,mask=mymask))
        self.assertTrue(len(rows),mymask.sum())
        for i,row in zip(mymask.nonzero()[0],rows):
            # make sure rows are correct length
            self.assertEquals(len(row),len(self.table1.keys()))
            
            # make sure cells are correct
            for j, col in enumerate(keyorder):
                if not isinstance(row[j],float) or (isinstance(row[j],float) and not numpy.isnan(row[j])):
                    val1 = row[j]
                    val2 = self.table1[col][i]
                    val3 = self.val_dict[col][i]
                    self.assertEquals(val1,val2,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val2))
                    self.assertEquals(val1,val3,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val3))

        # fetch masking columns ----------------------------
        some_keys = keyorder[:3]
        rows = list(self.table1.get_rows(keyorder=some_keys))
        self.assertTrue(len(rows),len(self.table1))
        for i,row in enumerate(rows):
            # make sure rows are correct length
            self.assertEquals(len(row),len(some_keys))
            
            # make sure cells are correct
            for j, col in enumerate(some_keys):
                if not isinstance(row[j],float) or (isinstance(row[j],float) and not numpy.isnan(row[j])):
                    val1 = row[j]
                    val2 = self.table1[col][i]
                    val3 = self.val_dict[col][i]
                    self.assertEquals(val1,val2,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val2))
                    self.assertEquals(val1,val3,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val3))

        # fetch masking columns and rows -------------------
        some_keys = keyorder[:3]
        rows = list(self.table1.get_rows(keyorder=some_keys,mask=mymask))
        self.assertTrue(len(rows),mymask.sum())
        for i,row in zip(mymask.nonzero()[0],rows):
            # make sure rows are correct length
            self.assertEquals(len(row),len(some_keys))
            
            # make sure cells are correct
            for j, col in enumerate(some_keys):
                if not isinstance(row[j],float) or (isinstance(row[j],float) and not numpy.isnan(row[j])):
                    val1 = row[j]
                    val2 = self.table1[col][i]
                    val3 = self.val_dict[col][i]
                    self.assertEquals(val1,val2,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val2))
                    self.assertEquals(val1,val3,
                                      "Row cell value does not match ArrayTable cell value (col %s): %s vs %s" % (col,val1,val3))

    def test_as_array(self):
        # pull only numerical keys, since numpy arrays can't be heterogeneous in dtype
        keyorder = [X[0] for X in self.table1.items() if X[1].dtype.kind in _NUMERIC_DTYPES]
        self.assertGreater(len(keyorder),0)

        # fetch whole table --------------------------------
        data = self.table1.as_array(keyorder=keyorder)
        self.assertTrue(data.shape[0],len(self.table1))
        self.assertTrue(data.shape[1],len(keyorder))
        
        for j,col in enumerate(keyorder):
            self.assertTrue(equal_enough(data[:,j],self.table1[col]))
            self.assertTrue(equal_enough(data[:,j],self.val_dict[col]))

        # fetch masking rows -------------------------------
        mymask  = self.table1["bool1"] == True
        mymask &= self.table1["bool2"] == True
        mymask |= self.table1["int1"] > 70000
        self.assertGreater(mymask.sum(),0)

        data = self.table1.as_array(keyorder=keyorder,mask=mymask)
        self.assertTrue(data.shape[0],mymask.sum())
        self.assertTrue(data.shape[1],len(keyorder))

        for j,col in enumerate(keyorder):
            self.assertTrue(equal_enough(data[:,j],self.table1[col][mymask]))
            self.assertTrue(equal_enough(data[:,j],self.val_dict[col][mymask]))

    def test_to_from_file(self):
        """Test read/write from files of float, int, bool, and str columns,
        with and without row and column selection
        """
        # write out, read in whole table --------------------------------
        fout = tempfile.NamedTemporaryFile("w",delete=False)
        self.table1.to_file(fout)
        fout.close()
        table1_reloaded = ArrayTable.from_file(open(fout.name))

        self.assertSetEqual(set(self.table1.keys()),set(table1_reloaded.keys()))
        for k in self.table1.keys():
            msg = "Failed read/write equality test for %s" % k
            self.assertTrue(equal_enough(self.table1[k],table1_reloaded[k],printer=printer),msg)

        os.remove(fout.name)

        # write out, read in, and test, masking rows -------------------
        mymask  = self.table1["bool1"] == True
        mymask &= self.table1["bool2"] == True
        mymask |= self.table1["int1"] > 70000
        self.assertGreater(mymask.sum(),0)
        
        fout = tempfile.NamedTemporaryFile("w",delete=False)
        self.table1.to_file(fout,mask=mymask)
        fout.close()
        table1_reloaded = ArrayTable.from_file(open(fout.name))
        
        # make sure rows are all there
        self.assertEquals(len(table1_reloaded),mymask.sum())
        
        # make sure columns are all there
        self.assertSetEqual(set(table1_reloaded.keys()),set(self.table1.keys()))
        
        # make sure values are correct
        for k in self.table1.keys():
            self.assertTrue(equal_enough(self.table1[k][mymask],table1_reloaded[k]))
        
        os.remove(fout.name)
        
        # write out, read in, and test, masking columns ----------------
        keys = ["string1","bool1","int2"]

        fout = tempfile.NamedTemporaryFile("w",delete=False)
        self.table1.to_file(fout,keyorder=keys)
        fout.close()
        table1_reloaded = ArrayTable.from_file(open(fout.name))
        
        # make sure rows are all there
        self.assertEquals(len(table1_reloaded),len(self.table1))

        # make sure columns are all there
        self.assertSetEqual(set(table1_reloaded.keys()),set(keys))
        
        # make sure values are correct
        for k in keys:
            msg = "Failure in recovery from file for key %s" % k
            self.assertTrue(equal_enough(self.table1[k],table1_reloaded[k]),msg)
        
        os.remove(fout.name)
        
        # write out, read in, and test, masking columns and rows -------
        fout = tempfile.NamedTemporaryFile("w",delete=False)
        self.table1.to_file(fout,keyorder=keys,mask=mymask)
        fout.close()
        table1_reloaded = ArrayTable.from_file(open(fout.name))
        
        # make sure rows are all there
        self.assertEquals(len(table1_reloaded),mymask.sum())

        # make sure columns are all there
        self.assertSetEqual(set(table1_reloaded.keys()),set(keys))
        
        # make sure values are correct
        for k in keys:
            self.assertTrue(equal_enough(self.table1[k][mymask],table1_reloaded[k]))
        
        os.remove(fout.name)
    
    def test_as_data_frame(self):
        t2 = copy.deepcopy(self.table1)
        df = t2.as_data_frame()
        
        # check identity of keys
        self.assertSetEqual(set(df.columns),set(self.table1.keys()))
        
        # check identity of values
        for k in self.table1.keys():
            equal_enough(self.table1[k],df[k])
        
    def test_from_file_to_data_frame(self):
        fout = tempfile.NamedTemporaryFile("w",delete=False)
        self.table1.to_file(fout)
        fout.close()
        df = ArrayTable.from_file_to_data_frame(open(fout.name))
        
        # check identity of keys
        self.assertSetEqual(set(df.columns),set(self.table1.keys()))
        
        # check identity of values
        for k in self.table1.keys():
            equal_enough(self.table1[k],df[k])
        
        os.remove(fout.name)

    def test_concat(self):
        new_table = ArrayTable.concat(self.table1,self.table1)
        
        # check length
        self.assertEquals(len(new_table),2*len(self.table1))
        
        # check columns
        self.assertSetEqual(set(new_table.keys()),set(self.table1.keys()))
        
        # check values
        for k in new_table:
            self.assertTrue(equal_enough(new_table[k][:len(self.table1)],self.table1[k]))
            self.assertTrue(equal_enough(new_table[k][len(self.table1):],self.table1[k]))
        
        # make sure we error out if columns are different
        prob_table = copy.deepcopy(self.table1)
        prob_table["nope"] = numpy.ones(len(prob_table))
        self.assertRaises(ValueError,ArrayTable.concat,self.table1,prob_table)

    def test_merge(self):
        keys = list(self.table1.keys())
        keys1 = keys[:3]
        keys2 = keys[3:]
        table1 = ArrayTable({ K : copy.deepcopy(self.val_dict[K]) for K in keys1 })
        table2 = ArrayTable({ K : copy.deepcopy(self.val_dict[K]) for K in keys2 })
        
        table1["merge_key"] = numpy.arange(0,len(table1))
        table2["merge_key"] = numpy.arange(0,len(table2))
        
        table1.merge(table2,"merge_key","merge_key")
        self.assertSetEqual(set(table1.keys()) - {"merge_key"},set(self.table1.keys()))
        for k in self.table1.keys():
            self.assertTrue(equal_enough(self.table1[k],table1[k]))

    def test_equal_enough(self):
        # positive & negative controls for floats
        for i in range(-10,10):
            a = 10**i * numpy.random.randn(100000)
            b = copy.deepcopy(a)
            if (a != 0).any():
                self.assertTrue(equal_enough(a,a))
                self.assertTrue(equal_enough(a,b))
                numpy.random.shuffle(b)
                self.assertFalse(equal_enough(a,b))
        
        # positive & negative controls for ints
        a = numpy.arange(-100000,100000,5)
        b = copy.deepcopy(a)
        self.assertTrue(equal_enough(a,a))
        self.assertTrue(equal_enough(a,a.astype(float)))
        self.assertTrue(equal_enough(a,b))
        numpy.random.shuffle(b)
        self.assertFalse(equal_enough(a,b))
        
        # positive & negative controls for floats at various levels of tolerance
        tol = 1e-10
        for my_root in numpy.linspace(-30,2,5):
            for my_tol in numpy.arange(-16,0):
                a = (10**my_root)*numpy.random.random(size=1000)
                noise = (10**my_tol)*numpy.random.random(size=1000)
                if (noise > tol).any():
                    msg = "Failed expected false tolerance test at %s noise level %s" % (my_root,10**my_tol)
                    self.assertFalse(equal_enough(a,a+noise,tol=tol),msg)
                else:
                    msg = "Failed expected true tolerance test at %s noise level %s" % (my_root,10**my_tol)
                    self.assertTrue(equal_enough(a,a+noise,tol=tol),msg)
        
        # positive and negative controls for nans
        a = numpy.random.random(size=100000)
        b = copy.deepcopy(a)
        self.assertTrue((a == b).all())
        self.assertTrue(equal_enough(a,b))
        indices = numpy.random.randint(0,high=len(a),size=50)
        a[indices] = numpy.nan
        self.assertTrue(equal_enough(a,a))
        self.assertFalse(equal_enough(a,b))
        b[indices] = numpy.nan
        self.assertTrue(equal_enough(a,b))
        
        # positive and negative controls for infs
        a = numpy.random.random(size=100000)
        b = copy.deepcopy(a)
        self.assertTrue((a == b).all())
        self.assertTrue(equal_enough(a,b))
        
        indices1 = numpy.random.randint(0,high=len(a),size=50)
        indices2 = numpy.random.randint(0,high=len(a),size=50)
        
        a[indices1] = numpy.inf
        self.assertTrue(equal_enough(a,a))
        self.assertFalse(equal_enough(a,b))
        
        b[indices1] = -numpy.inf
        self.assertFalse(equal_enough(a,b))
        
        b[indices1] = numpy.inf
        self.assertTrue(equal_enough(a,b))
        
        a[indices2] = -numpy.inf
        self.assertFalse(equal_enough(a,b))
        
        b[indices2] = -numpy.inf
        self.assertTrue(equal_enough(a,b))
        
        # positive & negative controls for dtype test
        a = numpy.random.random(size=100000)
        b = numpy.array([str(X) for X in a])
        self.assertFalse(equal_enough(a,b))
        self.assertTrue(equal_enough(a,a))
        self.assertTrue(equal_enough(b,b))

        
