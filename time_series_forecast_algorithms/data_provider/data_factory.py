from data_provider.data_loader import Dataset_Custom, Dataset_Custom_2, Dataset_Pred, Dataset_PEMS, Dataset_solar, Dataset_MGV
from torch.utils.data import DataLoader

data_dict = {
    'custom': Dataset_Custom,
    'AEP': Dataset_Custom,
    'air_pollution': Dataset_Custom,
    'AirQuality': Dataset_Custom,
    'AUROPHARMA': Dataset_Custom,
    'AXISBANK': Dataset_Custom,
    'BSD': Dataset_Custom,
    'BTCBUSD': Dataset_Custom,
    'BTCUSDT': Dataset_Custom,
    'Chicago_weather': Dataset_Custom,
    'daily_sunspots': Dataset_Custom,
    'DST': Dataset_Custom,
    'ECG': Dataset_Custom,
    'electricity': Dataset_Custom,
    'electricity_usage': Dataset_Custom,
    'ETTh1': Dataset_Custom,
    'ETTh2': Dataset_Custom,
    'ETTm1': Dataset_Custom,
    'ETTm2': Dataset_Custom,
    'exchange_rate': Dataset_Custom,
    'GEFCom2014-E': Dataset_Custom_2,
    'GEFCom2014-P': Dataset_Custom,
    'GEFCom2014-S': Dataset_Custom,
    'global_temp': Dataset_PEMS,
    'global_wind': Dataset_PEMS,
    'gold_rate': Dataset_Custom,
    'grid_loss': Dataset_Custom,
    'HHNGSP': Dataset_Custom,
    'HPC': Dataset_Custom,
    'JanataHack': Dataset_Custom,
    'MGV_eccentricity': Dataset_MGV,
    'MGV_missing_tooth': Dataset_MGV,
    'MGV_no_fault': Dataset_MGV,
    'MITV': Dataset_Custom,
    'MPW': Dataset_Custom,
    'MLTemp': Dataset_Custom,
    'NEL': Dataset_Custom,
    'NGPF': Dataset_Custom,
    'PEMS03': Dataset_PEMS,
    'PEMS04': Dataset_PEMS,
    'PEMS07': Dataset_PEMS,
    'PEMS08': Dataset_PEMS,
    'solar_AL': Dataset_solar,
    'TCPC': Dataset_Custom,
    'temp_humi': Dataset_Custom,
    'traffic': Dataset_Custom,
    'TSG_IT': Dataset_Custom,
    'weather': Dataset_Custom,
    'WTH': Dataset_Custom,
    'WTPG1': Dataset_Custom_2,
    'WTPG2': Dataset_Custom_2,
    'WTPG3': Dataset_Custom_2,
    'WTPG4': Dataset_Custom_2
}

def data_provider(args, flag):
    Data = data_dict[args.data]
    timeenc = 0 if args.embed != 'timeF' else 1

    if flag == 'test':
        shuffle_flag = False
        drop_last = False  # fix bug
        batch_size = args.batch_size
        freq = args.freq
    elif flag == 'pred':
        shuffle_flag = False
        drop_last = False
        batch_size = 1
        freq = args.freq
        Data = Dataset_Pred
    else:
        shuffle_flag = True
        drop_last = True
        batch_size = args.batch_size
        freq = args.freq

    data_set = Data(
        root_path=args.root_path,
        data_path=args.data_path,
        flag=flag,
        size=[args.seq_len, args.label_len, args.pred_len],
        features=args.features,
        target=args.target,
        timeenc=timeenc,
        freq=freq,
        enc_in=args.enc_in
    )

    print(flag, len(data_set))
    data_loader = DataLoader(
        data_set,
        batch_size=batch_size,
        shuffle=shuffle_flag,
        num_workers=args.num_workers,
        drop_last=drop_last)
    return data_set, data_loader
