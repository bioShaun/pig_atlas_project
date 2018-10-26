from __future__ import print_function
import sys

if not len(sys.argv) == 4:
    msg = '    python ' + sys.argv[0] + ' table col.list reorder.table'
    print(msg)
    sys.exit(0)

input_table = sys.argv[1]
col_list_file = sys.argv[2]
output_table = sys.argv[3]

col_list = [each.strip() for each in open(col_list_file)]
out_col_list = []
output_table_inf = open(output_table, 'w')
with open(input_table) as input_table_inf:
    index_list = []
    for n, eachline in enumerate(input_table_inf):
        eachline_inf = eachline.rstrip().split('\t')
        if n == 0:
            for each in col_list:
                if each in eachline_inf:
                    index_list.append(eachline_inf.index(each))
                    out_col_list.append(each)
            output_table_inf.write('%s\t%s\n' % (
                eachline_inf[0], '\t'.join(out_col_list)))
        else:
            output_list = [eachline_inf[0]]
            for each in index_list:
                output_list.append(eachline_inf[each])
            output_table_inf.write('%s\n' % '\t'.join(output_list))
output_table_inf.close()
