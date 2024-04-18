import re
import gong_shtc


# pattern from
# https://nispdata.nso.edu/webProdDesc2/presenter.php?file=halpha_gong_filenames_overview.html
name_pattern = re.compile(
    """
    (?P<obs_site>bb|ct|le|ml|mr|td|ud)  # observing site
    (?P<prod_code>\w{3})                # product code (should be bqc)
    (?P<year>\d{2})                     # 2 digit year
    (?P<month>\d{2})                    # month
    (?P<day>\d{2})                      # date
    t                                   # literal
    (?P<hour>\d{2})                     # hour
    (?P<minute>\d{2})                   # minute
    (?P<second>\d{2})?                  # second, optional
    c                                   # literal, begins carrington rotation
    (?P<car_rot>\d{4})                  # carrington rotation number
    _                                   # literal, begins carrington longitude
    (?P<car_lot>\d{3})                  # carrington longitude at left edge

    (?:                                 # non-captured group, optional
    _(?P<flavor_tags>\w+)               # flavor tags
    )?
    """,
    re.MULTILINE)

# get current filename
# get list of filenames in series
# sort filenames
# pick closest one in time
