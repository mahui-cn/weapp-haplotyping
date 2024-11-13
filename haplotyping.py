# -*- coding: utf-8 -*-
import os
import re
import json
import logging

# logging.basicConfig(level=logging.INFO)


# Y/mtDNA单倍群分型
class Haplotyping:
    # 单倍群分型结果
    __haplogroup_list = None
    # 单倍群分型树
    __haplo_tree = None
    # 单倍群分型树的时间戳
    __timestamp = None
    # 单倍群分型树的来源
    __source = ""
    # 单倍群分型树是Y或mt
    __is_y_mt = ""
    # 单倍群分型树的单倍群总数
    __total_haplo_count = 0
    # 单倍群分型树的SNP总数
    __total_snp_count = 0
    # 单倍群分型路径中最少要确认的有连续derived SNP的单倍群数量，如果分型有此连续数量的单倍群有derived SNP，则认为是确定的分型结果
    __confirmed_positive_haplo = 0
    # 单倍群分型上游最多允许的没有derived SNP的单倍群数量，-1表示允许任何假阳SNP。如果分型后上游有超过此连续数量的单倍群没有derived SNP，则认为分型结果是假阳
    __allowed_negative_haplo = 0
    # 输出的单倍群分型数量
    __max_haplo_count = 0
    # 单倍群分型树的haplo键名
    __haplo_key = ""
    # 单倍群分型树的children键名
    __children_key = ""
    # 单倍群分型树的snp列表键名
    __snp_list_key = ""
    # 单倍群分型树的snp键名
    __snp_key = ""
    # 单倍群分型树的snp pos19键名
    __pos19_key = ""
    # 单倍群分型树的snp pos38键名
    __pos38_key = ""
    # 单倍群分型树的snp pos键名
    __pos_key = ""
    # 单倍群分型树的ancestral突变键名
    __ancestral_key = ""
    # 单倍群分型树的derived突变键名
    __derived_key = ""
    # 单倍群分型树的用户突变键名
    __user_geno_key = "u"

    @property
    def HaploTree(self):
        return self.__haplo_tree

    @property
    def Source(self):
        return self.__source

    @property
    def IsYorMt(self):
        return self.__is_y_mt.upper()

    @property
    def Timestamp(self):
        return self.__timestamp

    @property
    def HaploCount(self):
        return self.__total_haplo_count

    @property
    def SNPCount(self):
        return self.__total_snp_count

    @property
    def ConfirmedPositiveHaplo(self):
        return self.__confirmed_positive_haplo

    @property
    def AllowedNegativeHaplo(self):
        return self.__allowed_negative_haplo

    @property
    def MaxHaploCount(self):
        return self.__max_haplo_count

    @property
    def HaplogroupList(self):
        return self.__haplogroup_list

    def __init__(
        self,
        haploTreeFileName=None,
        source=None,
        isYorMt=None,
        confirmedPositiveHaplo=3,
        allowedNegativeHaplo=2,
        maxHaploCount=5,
        haploKey="n",
        childrenKey="c",
        snpListKey="m",
        snpKey="v",
        pos19Key="p19",
        pos38Key="p38",
        posKey="p",
        ancestralKey="a",
        derivedKey="d",
    ):
        self.__source = source
        self.__is_y_mt = isYorMt
        self.__confirmed_positive_haplo = confirmedPositiveHaplo
        self.__allowed_negative_haplo = allowedNegativeHaplo
        self.__max_haplo_count = maxHaploCount
        self.__haplo_key = haploKey
        self.__children_key = childrenKey
        self.__snp_list_key = snpListKey
        self.__snp_key = snpKey
        self.__pos19_key = pos19Key
        self.__pos38_key = pos38Key
        self.__pos_key = posKey
        self.__ancestral_key = ancestralKey
        self.__derived_key = derivedKey

        if haploTreeFileName == None:
            raise Exception("请指定单倍群树文件名")

        if not os.access(haploTreeFileName, os.F_OK):
            raise Exception("单倍群树文件不可访问：" + haploTreeFileName)

        if source == None:
            raise Exception("请指定单倍群树数据源")

        if not re.match("y", isYorMt, re.IGNORECASE) and not re.match(
            "mt", isYorMt, re.IGNORECASE
        ):
            raise Exception("请指定单倍群树是：y或mt")

        # 加载单倍群分型树
        with open(haploTreeFileName, "r", encoding="utf-8-sig") as haplo_tree_file:
            haplo_tree_json = json.load(haplo_tree_file)
            # 单倍群树的时间戳
            if "timestamp" in haplo_tree_json:
                self.__timestamp = haplo_tree_json["timestamp"]
            # 单倍群分型树是否在tree属性
            self.__haplo_tree = (
                haplo_tree_json["tree"]
                if "tree" in haplo_tree_json
                else haplo_tree_json
            )
        if self.__haplo_tree == None:
            raise Exception("单倍群树文件为空：" + haploTreeFileName)

    def __del__(self):
        self.__haplo_tree = None
        self.__haplogroup_list = None

    # 单倍群分析
    def analyse(self, user_genome: dict, genome_ref: str = "hg19") -> list:
        if user_genome == None or len(user_genome) == 0:
            raise Exception("用户基因数据为空")

        if not re.match("hg19", genome_ref, re.IGNORECASE) and not re.match(
            "hg38", genome_ref, re.IGNORECASE
        ):
            raise Exception("请指定用户基因数据的参考基因组是：hg19或hg38")

        # 遍历单倍群分型树
        self.__haplogroup_list = []
        self.__check_snp(self.__haplo_tree, user_genome, genome_ref, [])

        return self.__haplogroup_list

    # 递归单倍群树，检测用户每个SNP的突变情况。所有参数必须显式赋值，不能使用参数默认值，否则在多线程中，参数的默认值会在进程中共享，导致数据混乱
    def __check_snp(
        self,
        tree_node: dict,
        user_genome: dict,
        genome_ref: str,
        root_end_node_list: list,
    ):
        # 累计单倍群数量
        self.__total_haplo_count += 1

        # 记录从终端节点到根节点的路径，每个节点只记录一次
        if tree_node[self.__haplo_key] not in [
            haplo_node[self.__haplo_key] for haplo_node in root_end_node_list
        ]:
            root_end_node_list.insert(0, tree_node)

        if self.__snp_list_key in tree_node and len(tree_node[self.__snp_list_key]) > 0:
            # 累计SNP数量
            self.__total_snp_count += len(tree_node[self.__snp_list_key])

            # 检测每个SNP的突变情况
            for snp_dict in tree_node[self.__snp_list_key]:
                pos = ""
                if re.match("y", self.__is_y_mt, re.IGNORECASE):
                    if re.match("hg19", genome_ref, re.IGNORECASE):
                        pos = str(snp_dict[self.__pos19_key])
                    elif re.match("hg38", genome_ref, re.IGNORECASE):
                        pos = str(snp_dict[self.__pos38_key])
                elif re.match("mt", self.__is_y_mt, re.IGNORECASE):
                    pos = str(snp_dict[self.__pos_key])

                # 如果用户检测了此SNP，把用户突变值放在树上
                if pos in user_genome:
                    snp_dict[self.__user_geno_key] = user_genome[pos][0]
                    if snp_dict[self.__user_geno_key] == snp_dict[self.__derived_key]:
                        logging.info(
                            "\t{} 单倍群 {}，SNP位点 {} 产生突变：{} -> {}".format(
                                self.__is_y_mt.upper(),
                                tree_node[self.__haplo_key],
                                snp_dict[self.__snp_key],
                                snp_dict[self.__ancestral_key],
                                snp_dict[self.__derived_key],
                            )
                        )
                elif self.__user_geno_key in snp_dict:
                    # 如果用户未检测此SNP，但还有其他用户SNP值，则删除，避免数据混淆影响
                    del snp_dict[self.__user_geno_key]

        # 递归子节点
        if self.__children_key in tree_node and len(tree_node[self.__children_key]) > 0:
            for child_node in tree_node[self.__children_key]:
                self.__check_snp(
                    child_node,
                    user_genome,
                    genome_ref,
                    root_end_node_list,
                )

        else:  # 此节点没有子节点，是终端节点，开始计算单倍群分型

            # 当前单倍群路径上的每个单倍群和突变情况
            haplo_path_list = []
            # 此单倍群路径的单倍群分型结果
            haplogroup = None
            # 单倍群分型深度
            haplo_depth = 0
            # 单倍群路径中的累计用户derived SNP突变数
            total_user_der_count = 0
            # 单倍群路径中的累计用户SNP位点数
            total_user_var_count = 0
            # 测试覆盖的单倍群节点数
            tested_haplo = 0
            # 单倍群分型后出现的存在连续derived SNP单倍群数
            positive_haplo_count = 0
            # 单倍群分型后出现的存在连续无derived SNP单倍群数
            negative_haplo_count = 0

            # 遍历从终端节点到根节点的每个单倍群
            for haplo_node in root_end_node_list:
                if self.__snp_list_key in haplo_node:
                    # 当前单倍群节点中的用户已检测SNP数
                    user_var_count = 0
                    # 当前单倍群节点中的用户突变SNP数
                    user_der_count = 0

                    # 遍历当前单倍群节点的SNP突变情况，作为分型评分依据
                    for snp_dict in haplo_node[self.__snp_list_key]:
                        if self.__user_geno_key in snp_dict:
                            # 用户已检测SNP数
                            user_var_count += 1
                            if (
                                snp_dict[self.__user_geno_key]
                                == snp_dict[self.__derived_key]
                            ):
                                # 用户突变SNP数
                                user_der_count += 1

                    # 如果用户在当前单倍群节点有derived突变
                    if user_der_count > 0:
                        # 如果还未分型，则以此单倍群暂定分型（后续可能会因为上游节点存在无突变的单倍群节点，按照分型规则判断此处是假阳，则取消此分型结果）
                        if haplogroup == None:
                            haplogroup = haplo_node[self.__haplo_key]

                        # 连续阳性单倍群数量+1
                        positive_haplo_count += 1

                        # 连续阴性单倍群数量重置0
                        negative_haplo_count = 0

                        # 累计用户SNP位点汇总数
                        total_user_var_count += user_var_count

                        # 累计用户所有单倍群的derived SNP突变数
                        total_user_der_count += user_der_count

                    else:  # 如果用户在当前单倍群节点没有derived突变
                        # 判断下游连续阳性突变单倍群数量是否小于阈值
                        if positive_haplo_count < self.__confirmed_positive_haplo:
                            # 连续阳性单倍群数量重置0
                            positive_haplo_count = 0

                            # 连续阴性单倍群数量+1
                            negative_haplo_count += 1

                            # 如果当前单倍群节点的下游连续阴性单倍群数量大于阈值，且已经分型，则判断分型结果是跳变假阳，取消此前的分型结果
                            if (
                                self.__allowed_negative_haplo != -1
                                and negative_haplo_count > self.__allowed_negative_haplo
                                and haplogroup != None
                            ):
                                haplogroup = None
                                haplo_depth = 0
                                tested_haplo = 0
                                total_user_var_count = 0
                                total_user_der_count = 0

                    # 如果此时已分型
                    if haplogroup != None:
                        # 累计分型深度
                        haplo_depth += 1

                        # 如果当前单倍群节点中有用户已检测的SNP，则累计已检测的单倍群数
                        if user_var_count > 0:
                            tested_haplo += 1

                    # 保存此节点的单倍群名和SNP列表
                    haplo_path_list.append(
                        {
                            "haplo": haplo_node[self.__haplo_key],
                            "mutation": haplo_node[self.__snp_list_key],
                        }
                    )

            # 如果此单倍群分型路径有用户derived突变，且此单倍群分型结果在结果集中不存在，则新增
            if (
                total_user_der_count > 0
                and len(
                    [
                        haploObj
                        for haploObj in self.__haplogroup_list
                        if haploObj["haplo"] == haplogroup
                    ]
                )
                == 0
            ):
                # 根据此单倍群分型路径突变情况，计算分型结果可靠性评分
                haplo_score = 0
                if total_user_var_count != 0 and haplo_depth != 0:
                    # 用户derived SNP突变数/用户所有检测SNP数 * 有突变的单倍群节点数/单倍群分型深度
                    haplo_score = (total_user_der_count / total_user_var_count) * (
                        tested_haplo / haplo_depth
                    )

                # 加入单倍群分型结果列表
                self.__haplogroup_list.append(
                    {
                        "haplo": haplogroup,
                        "snp_derived_count": total_user_der_count,
                        "haplo_depth": haplo_depth,
                        "haplo_score": haplo_score,
                        "haplo_path": haplo_path_list,
                    }
                )

                # 按照规则排序单倍群分型结果
                self.__haplogroup_list.sort(
                    key=lambda haplo: (
                        haplo["snp_derived_count"],
                        haplo["haplo_depth"],
                        haplo["haplo_score"],
                    ),
                    reverse=True,
                )

                # 只保留指定数量的分型结果
                if len(self.__haplogroup_list) > self.__max_haplo_count:
                    self.__haplogroup_list.pop()

        # 处理完成每个节点后，删除这个节点，回到上层递归后，再压入下一个节点（兄弟节点或子节点）
        root_end_node_list.pop(0)

    # 输出单倍群分型结果，HTML表格
    def __str__(self):
        if self.__haplogroup_list == None or len(self.__haplogroup_list) == 0:
            return ""

        haploTable = []
        haploTable.append("<div class='table-responsive'>")
        haploTable.append("<table class='table'>")
        haploTable.append(
            "<thead><tr><th>{} 单倍群</th><th>SNP突变数</th><th>层级</th><th>可信度</th></tr></thead>".format(
                self.__is_y_mt.upper()
            )
        )
        haploTable.append("<tbody>")
        for idx, haplo in enumerate(self.__haplogroup_list):
            haploTable.append(
                "<tr style='{}'><td><a href='https://geneu.xyz/haplo-tree/{}/{}/{}' target='_blank'>{}</a></td><td>{}</td><td>{}</td><td>{:.2%}</td></tr>".format(
                    "color: red; font-size: larger;" if idx == 0 else "",
                    self.__source,
                    self.__is_y_mt,
                    haplo["haplo"],
                    (
                        "<span style='color: red; font-size: larger;'>{}</span>".format(
                            haplo["haplo"]
                        )
                        if idx == 0
                        else "<span style='color: white;'>{}</span>".format(
                            haplo["haplo"]
                        )
                    ),
                    haplo["snp_derived_count"],
                    haplo["haplo_depth"],
                    haplo["haplo_score"],
                )
            )
        haploTable.append("</tbody>")
        haploTable.append("</table>")
        haploTable.append("</div>")
        return "".join(haploTable)
