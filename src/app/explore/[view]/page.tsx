import { notFound } from "next/navigation";
import ExploreClient from "./explore-client";
import type { ViewKind, ClusterMeta } from "@/lib/types";
import proteinClusters from "@/../public/atlas/protein_clusters.json";
import moleculeClusters from "@/../public/atlas/molecule_clusters.json";
import natsynClusters from "@/../public/atlas/natsyn_clusters.json";

const VIEWS: ViewKind[] = ["protein", "molecule", "natsyn"];

const CLUSTERS_FOR: Record<ViewKind, ClusterMeta[]> = {
  protein: proteinClusters as ClusterMeta[],
  molecule: moleculeClusters as ClusterMeta[],
  natsyn: natsynClusters as ClusterMeta[],
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
  return <ExploreClient kind={kind} clusters={CLUSTERS_FOR[kind]} />;
}
