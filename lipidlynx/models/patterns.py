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

# lynx FA
fa_rgx = re.compile(r"(?P<link>[OP]-|FA|SBP)?(?P<c>\d{1,2})(:)(?P<db>\d)(?P<mod>[<].*[>])?")

# lynx PL
pl_rgx = re.compile(r"(?P<subclass>L)?(?P<class>P[ACEGIS]|PIP[1-3]?)(\s*\(?)?"
                    r"(?P<fa1>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position1>[_/\\])?"
                    r"(?P<fa2>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(\))")

# lynx SP
sp_rgx = re.compile(r"(?P<subclass>L)?(?P<class>Cer|SM|SP)(\s*\(?)?"
                    r"(?P<fa1>([OP]-|FA|SBP)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position1>[_/\\])?"
                    r"(?P<fa2>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(\))")

# lynx GL
gl_rgx = re.compile(r"(?P<fa1>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position1>[_/\\])?"
                    r"(?P<fa2>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position2>[_/\\])?"
                    r"(?P<fa3>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(\))")

# lynx CL
cl_rgx = re.compile(r"(?P<fa1>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position1>[_/\\])?"
                    r"(?P<fa2>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position2>[_/\\])?"
                    r"(?P<fa3>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(?P<position3>[_/\\])?"
                    r"(?P<fa4>([OP]-|FA)?(\d{1,2})(:)(\d)([<].*[>])?)?(\))")
