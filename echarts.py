# -*- coding: utf-8 -*-
import pyecharts.options as opts
from pyecharts.charts import Pie


# 使用祖源模型数据生成pie
def make_pie(admix_model):
    # 生成pie的坐标数据
    x_data = []
    y_data = []
    for admix in admix_model["admix"]:
        if round(admix["ratio"] * 100, 2) > 0:
            x_data.append(admix["name_cn"])
            y_data.append(admix["ratio"])

    data_pair = [list(z) for z in zip(x_data, y_data)]
    data_pair.sort(key=lambda x: x[1])

    (
        Pie(init_opts=opts.InitOpts(bg_color="#4fb1f7"))
        .add(
            series_name="祖源成分",
            data_pair=data_pair,
            rosetype="radius",
            radius="55%",
            center=["50%", "50%"],
            label_opts=opts.LabelOpts(is_show=False, position="center"),
        )
        .set_global_opts(
            title_opts=opts.TitleOpts(
                title=admix_model["name_cn"],
                subtitle=admix_model["desc_cn"],
                pos_left="center",
                pos_top="20",
                padding=5,
                title_textstyle_opts=opts.TextStyleOpts(color="#fff"),
            ),
            legend_opts=opts.LegendOpts(is_show=True),
        )
        .set_series_opts(
            tooltip_opts=opts.TooltipOpts(trigger="item", formatter="{b}: {d}%"),
            label_opts=opts.LabelOpts(color="rgba(255, 255, 255, 0.3)"),
        )
        .render("customized_pie.html")
    )
