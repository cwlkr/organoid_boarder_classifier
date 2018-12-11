import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse


def split(u, v, points):
    # return points on left side of UV
    return [p for p in points if np.cross(p - u, v - u) < 0]


def extend(u, v, points):
    if not points:
        return []

    # find furthest point W, and split search to WV, UW
    w = min(points, key=lambda p: np.cross(p - u, v - u))
    p1, p2 = split(w, v, points), split(u, w, points)
    return extend(w, v, p1) + [w] + extend(u, w, p2)


def convex_hull(points):
    # find two hull points, U, V, and split to left and right search
    u = min(points, key=lambda p: p[0])
    v = max(points, key=lambda p: p[0])
    left, right = split(u, v, points), split(v, u, points)

    # find convex hull on each side, last [v is duplicate!]
    return [v] + extend(u, v, left) + [u] + extend(v, u, right)


def radial():
    # radial classfier
    df = pd.read_csv('Nuclei.csv')
    xy = df[df.ImageNumber == 1][[
        'Location_Center_X', 'Location_Center_Y']].values
    hull = np.array(convex_hull(xy))

    x = np.array((hull[:, 0], hull[:, 1]))
    z = np.array((np.mean(xy[:, 0]), np.mean(xy[:, 1])))
    r = np.sqrt(np.sum((z[:, np.newaxis] - x)**2, axis=0))

    plt.plot(xy[:, 0], xy[:, 1], 'o')
    plt.plot(hull[:, 0], hull[:, 1], 'o', color='red')
    plt.plot(hull[:, 0], hull[:, 1], color='red')
    for t in xy:
        x = np.sqrt(np.sum((t-z)**2))
        if x > (np.min(r)*0.9):
            plt.plot(t[0], t[1], 'o', color='red')
    plt.show()


def boarder_classifier(points, thresh, imgnr, plot):
    # how to keep track of indexing

    hull = np.array(convex_hull(points))
    not_hull = np.array([np.array(p) for p in points if p not in hull])

    pair_dist = np.sqrt(
        np.sum(np.power(hull[:, :] - not_hull[:, np.newaxis, :], 2), axis=2))
    argsort = np.argsort(pair_dist, axis=1)
    AB_index = argsort[:, :2]
    AB = hull[AB_index, :]
    ab = np.sqrt(
        np.sum(np.power(AB[:, :] - not_hull[:, np.newaxis, :], 2), axis=2))
    a = ab[:, 0].flatten()
    b = ab[:, 1].flatten()
    c = np.sqrt(np.sum(np.diff(AB, axis=1)**2, axis=2)).flatten()
    cosa = (b**2 + c**2 - a**2) / (2 * b * c)
    h = b * np.sin(np.arccos(cosa))
    h < thresh
    if plot:
        plot_img(hull, not_hull, thresh, h, imgnr)
    # all should have same index!
    res = np.zeros(len(points))

    n_ind = 0
    res[np.all(np.isin(points, hull), axis=1)] = 1
    for i in range(len(points)):
        if np.all(np.isin(points, not_hull), axis=1)[i]:
            res[i] = 1 if h[n_ind] < thresh else 0
            n_ind = n_ind + 1
    return res


def plot_img(hull, not_hull, thresh, h, imgnr):
    fig, ax = plt.subplots()
    ax.plot(hull[:, 0], hull[:, 1], 'o', color='red')
    ax.plot(hull[:, 0], hull[:, 1], color='red')
    for i in range(len(not_hull)):
        if h[i] < thresh:
            ax.plot(not_hull[i, 0], not_hull[i, 1], 'o', color='red')
        else:
            ax.plot(not_hull[i, 0], not_hull[i, 1], 'o', color='blue')
    fig.savefig('fig' + str(imgnr) + '.pdf')
    plt.show()


if __name__ == '__main__':
    # classifier
    # all hull points are classified as boarder points
    # calc dist from boundery for each point not in boundery
    # for each point in not_hull calc dist to hull points
    desc = "Classifies the cells of a 2D slice of an Organiod according if they are a boarder cell or not"
    parser = argparse.ArgumentParser(description=desc)
    defaults = {'xlabel': 'Location_Center_X',
                'ylabel': 'Location_Center_Y',
                'arealabel': 'AreaShape_Area',
                'imagenr': 'ImageNumber',
                'colname': 'boarder_cell'}

    parser.add_argument(
        'input', type=str,
        help="absolute path to csv file, or filename if its in the same folder as the .py file.")
    parser.add_argument(
        'output', type=str,
        help='output path of the filename in obsolute or filename, if filename it gets saved to folder of .py file'
    )
    parser.add_argument(
        '--surpress_plot', '-p', action='store_false',
        default=True,
        help="If set, does not create Figures of the Cell classification"
    )
    parser.add_argument(
        '--xlabel', '-x', type=str,
        default=defaults['xlabel'],
        help="Manually set x label of csv. Default is \"Location_Center_X\""
    )
    parser.add_argument(
        '--ylabel', '-y', type=str,
        default=defaults['ylabel'],
        help="Manually set x label of csv. Default is \"Location_Center_Y\""
    )
    parser.add_argument(
        '--arealabel', '-a', type=str,
        default=defaults['arealabel'],
        help="Manually set label the area of a cell is stored in the csv. Default is \"AreaShape_Area\""
    )
    parser.add_argument(
        '--imagenr', '-i', type=str,
        default=defaults['imagenr'],
        help="Manually set label the unique identifier for the image is stored in the csv, deafult is  \"ImageNumber\""
    )
    parser.add_argument(
        '--colname', '-c', type=str,
        default=defaults['colname'],
        help="Specify name of column the boarder classification is saved to, default is \"boarder_cell\""
    )
    args = parser.parse_args()
    print(args)

    df = pd.read_csv(args.input)
    res = []
    len(df[args.imagenr].values)
    for i in range(len(pd.unique(df[args.imagenr].values))):
        points = df[df[args.imagenr] == (i+1)][[args.xlabel, args.ylabel]].values
        min_dmi = np.min(df[df[args.imagenr] == (i+1)][[args.arealabel]].values/np.pi)
        res.append(pd.DataFrame(boarder_classifier(points, min_dmi/2, i+1, args.surpress_plot)))
    df[args.colname] = pd.concat(res).values
    df.to_csv(args.output)

    ############################
    # test
    # points = np.array([(0, 0), (0.5, 0.5), (1, 1), (0, 1), (1, 0), (0.5, 0.01), (0.75, 0.75), (0.5,0.09)])
    # points
    # any(boarder_classifier(points, 0.1, 15, True) == np.array([1., 0., 1., 1., 1., 1., 0., 1.]))
