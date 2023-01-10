Python 3.10.7 (tags/v3.10.7:6cc6b13, Sep  5 2022, 14:08:36) [MSC v.1933 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.

====== RESTART: G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py =====
Structure dict
KeyboardInterrupt
parsecompositiontomatrix0()
enter 'lite' for dfb while others will be dfa
Traceback (most recent call last):
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\indexes\base.py", line 3800, in get_loc
    return self._engine.get_loc(casted_key)
  File "pandas\_libs\index.pyx", line 138, in pandas._libs.index.IndexEngine.get_loc
  File "pandas\_libs\index.pyx", line 165, in pandas._libs.index.IndexEngine.get_loc
  File "pandas\_libs\hashtable_class_helper.pxi", line 5745, in pandas._libs.hashtable.PyObjectHashTable.get_item
  File "pandas\_libs\hashtable_class_helper.pxi", line 5753, in pandas._libs.hashtable.PyObjectHashTable.get_item
KeyError: 'Structure dict'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "<pyshell#0>", line 1, in <module>
    parsecompositiontomatrix0()
  File "G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py", line 32, in parsecompositiontomatrix0
    dictparser = tmp['Structure dict']
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\series.py", line 982, in __getitem__
    return self._get_value(key)
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\series.py", line 1092, in _get_value
    loc = self.index.get_loc(label)
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\indexes\base.py", line 3802, in get_loc
    raise KeyError(key) from err
KeyError: 'Structure dict'

====== RESTART: G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py =====
parsecompositiontomatrix0()
enter 'lite' for dfb while others will be dfa
Traceback (most recent call last):
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\indexes\base.py", line 3800, in get_loc
    return self._engine.get_loc(casted_key)
  File "pandas\_libs\index.pyx", line 138, in pandas._libs.index.IndexEngine.get_loc
  File "pandas\_libs\index.pyx", line 165, in pandas._libs.index.IndexEngine.get_loc
  File "pandas\_libs\hashtable_class_helper.pxi", line 5745, in pandas._libs.hashtable.PyObjectHashTable.get_item
  File "pandas\_libs\hashtable_class_helper.pxi", line 5753, in pandas._libs.hashtable.PyObjectHashTable.get_item
KeyError: 'Structure dict'

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "<pyshell#2>", line 1, in <module>
    parsecompositiontomatrix0()
  File "G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py", line 32, in parsecompositiontomatrix0
    dictparser = tmp['Structure dict']
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\series.py", line 982, in __getitem__
    return self._get_value(key)
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\series.py", line 1092, in _get_value
    loc = self.index.get_loc(label)
  File "C:\Users\hnstseng\AppData\Local\Programs\Python\Python310\lib\site-packages\pandas\core\indexes\base.py", line 3802, in get_loc
    raise KeyError(key) from err
KeyError: 'Structure dict'

====== RESTART: G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py =====
parsecompositiontomatrix0()
enter 'lite' for dfb while others will be dfa
MS1scan no                                                                 14819
in [H+]                                                              1813.946262
Theoretical mass                                                       1813.9435
Structure dict                    {'F': '1', 'H': '3', 'N': '4', 'S': 0, 'G': 0}
Structure                                                                 F1H3N4
MS2 Scan no                                                                14857
hitpeaklist                    [344.17, 374.181, 376.196, 432.222, 450.233, 4...
normalizedselectedintensity    [5.005929610596618, 0.1689084519950698, 0.9342...
Name: 0, dtype: object
Enter the file namerevised_zfNGbrain_20221129_1219_proc1
Finish exporting revised_zfNGbrain_20221129_1219_proc1.csv

====== RESTART: G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py =====
parsecompositiontomatrix0()
                                
enter 'lite' for dfb while others will be dfa
Hello this is line 1, and we'll get dict keys
<class 'dict'> {'F': 0, 'H': '5', 'N': '2', 'S': 0, 'G': 0, 'Na': 1}
{'F': 0, 'H': '5', 'N': '2', 'S': 0, 'G': 0, 'Na': 1}
['F', 'H', 'N', 'S', 'G', 'Na']
Enter the file namerevised_zfNGbrain_20221129_1219_proc1
Finish exporting revised_zfNGbrain_20221129_1219_proc1.csv

====== RESTART: G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py =====
parseintensitytomatrix0()
                                
Enter the file namerevised_zfNGbrain_20221129_1219_proc2
Finish exporting revised_zfNGbrain_20221129_1219_proc2.csv
>>> 
====== RESTART: G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py =====
>>> fixhitlistbugmethodA()
...                                 
modified? ['Unnamed: 0.1', '344.169', '344.17', '344.171', '374.176', '374.178', '374.179', '374.18', '374.181', '374.182', '374.183', '376.188', '376.189', '376.19', '376.192', '376.195', '376.196', '376.197', '376.198', '376.199', '406.197', '406.201', '406.203', '406.204', '406.205', '406.206', '406.207', '406.208', '406.209', '406.21', '432.221', '432.222', '432.223', '432.224', '432.226', '450.231', '450.232', '450.233', '450.234', '450.235', '464.245', '464.246', '464.247', '464.248', '464.249', '464.25', '464.251', '464.252', '464.254', '580.293', '580.294', '580.295', '580.296', '580.297', '580.298', '580.299', '793.391', '793.393', '793.394', '793.395', '793.396', '793.397', '793.398', '793.399', '793.401', '793.404', '793.406', '793.408', '793.409', '825.404', '825.416', '825.418', '825.419', '825.42', '825.421', '825.422', '825.423', '825.424', '825.425', '825.426', '825.427', 'F', 'G', 'H', 'MS1scan no', 'MS2 Scan no', 'N', 'Na', 'Questionable', 'S', 'Structure', 'Structure dict', 'Theoretical mass', 'Unnamed: 0', 'hitpeaklist', 'in [H+]', 'normalizedselectedintensity', 344.1547, 374.1653, 376.1966, 432.2071, 450.2334, 464.249, 793.3808, 825.4227, 450.2334, 344.1547, 374.1653, 450.2334, 450.2334, 825.4227, 376.1966, 432.2071, 374.1653, 464.249, 793.3808, 825.4227, 344.1547, 450.2334, 374.1653, 376.1966, 406.2072, 825.4227, 793.3808, 793.3808, 825.4227, 825.4227, 432.2071, 825.4227, 406.2072, 793.3808, 376.1966, 406.2072, 793.3808, 793.3808, 825.4227, 580.2964, 580.2964, 580.2964, 825.4227, 580.2964, 793.3808, 406.2072, 464.249, 580.2964, 432.2071, 464.249, 374.1653, 793.3808, 374.1653, 464.249, 406.2072, 793.3808, 580.2964, 432.2071, 464.249, 464.249, 825.4227, 406.2072, 464.249, 825.4227, 376.1966, 825.4227, 793.3808, 406.2072, 406.2072, 376.1966, 376.1966, 376.1966, 406.2072, 406.2072, 376.1966, 580.2964, 464.249, 793.3808, 374.1653, 793.3808]

Warning (from warnings module):
  File "G:\GlycoMSParser\src\parsecomposition_for_analyticandML.py", line 156
    tmpopt.set_axis(c, axis=1, inplace=True)#Warning
FutureWarning: DataFrame.set_axis 'inplace' keyword is deprecated and will be removed in a future version. Use `obj = obj.set_axis(..., copy=False)` instead
modified index Index(['Unnamed: 0.1',      '344.169',       '344.17',      '344.171',
            '374.176',      '374.178',      '374.179',       '374.18',
            '374.181',      '374.182',
       ...
             376.1966,       376.1966,       406.2072,       406.2072,
             376.1966,       580.2964,        464.249,       793.3808,
             374.1653,       793.3808],
      dtype='object', length=177)
Index([                 'MS1scan no',                     'in [H+]',
                  'Theoretical mass',                   'Structure',
                    'Structure dict',                 'MS2 Scan no',
                       'hitpeaklist', 'normalizedselectedintensity',
                                 'F',                           'H',
                                 'N',                           'S',
                                 'G',                          'Na',
                      'Questionable',                      344.1547,
                            374.1653,                      376.1966,
                            406.2072,                      432.2071,
                            450.2334,                       464.249,
                            580.2964,                      793.3808,
                            825.4227],
      dtype='object')
Enter the file namerevised_zfNGbrain_20221129_1219_proc3
Finish exporting revised_zfNGbrain_20221129_1219_proc3.csv
