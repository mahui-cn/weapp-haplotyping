# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import io
import base64


# 生成饼图
def make_pie(pie_x, pie_label, file_format="png"):
    plt.style.use("_mpl-gallery-nogrid")

    # make data
    colors = plt.get_cmap("Blues")(np.linspace(0.2, 0.7, len(pie_x)))

    # plot
    fig, ax = plt.subplots()
    ax.pie(
        x=pie_x,
        labels=pie_label,
        colors=colors,
        radius=3,
        center=(4, 4),
        labeldistance=0.2,
        autopct="%1.2f%%",
        wedgeprops={"linewidth": 1, "edgecolor": "white"},
        frame=True,
    )

    ax.set(xlim=(0, 8), xticks=np.arange(1, 8), ylim=(0, 8), yticks=np.arange(1, 8))

    # 返回图片的base64编码的字符串
    image_buf = io.BytesIO()
    plt.savefig(image_buf, format=file_format)
    image_buf.seek(0)
    image_b64_str = base64.b64encode(image_buf.read()).decode()

    return image_b64_str
