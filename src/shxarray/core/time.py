# This file is part of frommle2.
# frommle2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# frommle2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with Frommle; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# Author Roelof Rietbroek (r.rietbroek@utwente.nl), 2022

from datetime import datetime,timedelta

def decyear2dt(decyear):
    """Convert a decimal year to a datetime object"""
    year=int(decyear)
    if year == 0:
        return datetime.min
    
    jan1=datetime(year,1,1)
    return jan1+(decyear-year)*(datetime(year+1,1,1)-jan1)


def dt2decyear(dt):
    """Convert a datetime object to a decimal year"""
    year=dt.year

    jan1=datetime(year,1,1)
    
    jan1next=datetime(year+1,1,1)
    yrlen=(jan1next-jan1).total_seconds()
    return year+(dt-jan1).total_seconds()/yrlen

