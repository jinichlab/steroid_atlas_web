import { notFound } from "next/navigation";
import ExploreClient from "./explore-client";
import type { ViewKind, ClusterMeta } from "@/lib/types";
import proteinClusters from "@/../public/atlas/protein_clusters.json";
import moleculeClusters from "@/../public/atlas/molecule_clusters.json";
import natsynClusters from "@/../public/atlas/natsyn_clusters.json";
// High-dim k-means clusters ([] until analysis/scripts/19_highdim_kmeans.py
// has run) — they drive the clustering toggle.
import proteinClustersHd from "@/../public/atlas/protein_clusters_hd.json";
import moleculeClustersHd from "@/../public/atlas/molecule_clusters_hd.json";

const VIEWS: ViewKind[] = ["protein", "molecule", "natsyn"];

const CLUSTERS_FOR: Record<ViewKind, ClusterMeta[]> = {
  protein: proteinClusters as ClusterMeta[],
  molecule: moleculeClusters as ClusterMeta[],
  natsyn: natsynClusters as ClusterMeta[],
};

const CLUSTERS_HD_FOR: Record<ViewKind, ClusterMeta[]> = {
  protein: proteinClustersHd as ClusterMeta[],
  molecule: moleculeClustersHd as ClusterMeta[],
  natsyn: [],
};

export function generateStaticParams() {
  return VIEWS.map((view) => ({ view }));
}

export default async function ExplorePage({
  params,
}: {
  params: Promise<{ view: string }>;
}) {
  const { view } = await params;
  if (!VIEWS.includes(view as ViewKind)) return notFound();
  const kind = view as ViewKind;
  return (
    <ExploreClient
      kind={kind}
      clusters={CLUSTERS_FOR[kind]}
      clustersHd={CLUSTERS_HD_FOR[kind]}
    />
  );
}
