import re

# general
rgx_blank = re.compile(r"\W+")


# lynx modification
# mod_rgx = re.compile(r'([{](\d{1,2}[EZ]?)(,\d{1,2}[EZ]?)*[}])?'
#                      r'(,?\d?\w{1,6}([{](\d{1,2}[RS]?)(,\d{1,2}[RS]?)*[}]))*')
mod_rgx = re.compile(r'((\d?\w{1,6})?([{][0-9EZRS,]{1,256}[}])?)')

mod_db_rgx = re.compile(r'[{](\d{1,2}[EZ]?)(,\d{1,2}[EZ]?)*[}]')

mod_general_rgx = re.compile(r',?(?P<count>\d)?(?P<cv>\w{1,6})'
                             r'[{](?P<positions>(\d{1,2}[RS]?)(,\d{1,2}[RS]?)*)[}]')

mod_position_rgx = re.compile(r'(?P<position>\d{1,2})(?P<position_type>[EZRS])?')
