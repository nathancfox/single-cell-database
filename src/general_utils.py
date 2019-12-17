import datetime as dt

def get_timestamp(mode = 'both'):
    now = dt.datetime.now()
    date_stamp = now.strftime('%y%m%d')
    time_stamp = now.strftime('%H%M%S')
    if mode == 'both':
        return(date_stamp + '_' + time_stamp)
    elif mode == 'date':
        return(date_stamp)
    elif mode == 'time':
        return(time_stamp)
    else:
        raise ValueError('mode must be one of {\'both\', \'date\', \'time\'}!')

def main():
    pass

if __name__ == '__main__':
    main()