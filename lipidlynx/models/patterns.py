import re

# general
rgx_blank = re.compile(r"\W+")


# lynx modification

mod_rgx = re.compile(r'((,)?([+-]?\d{0,2}\w{1,6})?([{][0-9EZRS,]{1,256}[}])?)')

mod_db_rgx = re.compile(r'[{](,?\d{1,2}[EZ]?)*[}]')

mod_lv0_delta_rgx = re.compile(r'(,?[-+]\d{1,2})+$')

mod_lv1_elem_rgx = re.compile(r',?(?P<count>[-+]\d{0,2})(?P<cv>[CHNOPSDK][a]?)+')  # elem can be Na

mod_lv2_types_rgx = re.compile(r',?(?P<count>\d{0,2})(?P<cv>\w{1,6})')

mod_lv3_position_rgx = re.compile(r',?(?P<count>\d{0,2})(?P<cv>\w{1,6})'
                                  r'[{](?P<positions>(\d{1,2}[RS]?)(,\d{1,2}[RS]?)*)[}]')

mod_lv4_position_rgx = re.compile(r'(?P<position>\d{1,2})(?P<position_type>[EZRS])?')

mod_sub_lv1_db_rgx = re.compile(r'[{](,?\d{1,2})*[}]')

mod_sub_lv2_db_rgx = re.compile(r'[{](,?\d{1,2}[EZ])*[}]')

# lynx FA

fa_rgx = re.compile(r'(?P<link>FA|[OP]-)?(?P<c>\d{1,2}):(?P<db>\d)(?P<mod><.*>)?')
