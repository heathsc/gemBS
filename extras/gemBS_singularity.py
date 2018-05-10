#!/usr/bin/python
# Hacked from easy install entry point script 
# so that we can sanitize sys.path, removing
# any entry that is no in /usr
__requires__ = 'gemBS==2.1.0'
import sys
# Clean sys.path
p=[]
for x in sys.path[:]:
	print(x)
	if x[:4] == '/usr':
		p.append(x)
sys.path = p
# Now continue launch of gemBS
from pkg_resources import load_entry_point

if __name__ == '__main__':
    sys.exit(
        load_entry_point('gemBS==2.1.0', 'console_scripts', 'gemBS')()
    )
