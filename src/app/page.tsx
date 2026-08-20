import Link from "next/link";
import summary from "@/../public/atlas/summary.json";
import type { Summary } from "@/lib/types";

const S = summary as Summary;

export default function Home() {
  return (
    <main className="mx-auto max-w-5xl px-6 py-14">
      <div className="mb-3 flex items-center gap-3 text-sm font-medium text-slate-500">
        <span className="text-lg">🧭</span> STEROID ATLAS
      </div>
      <h1 className="mb-6 text-4xl font-bold tracking-tight text-slate-900 sm:text-5xl">
        Explore steroids and the proteins they interact with.
      </h1>
      <p className="mb-8 max-w-2xl text-lg text-slate-600">
        An interactive UMAP atlas of steroid- and bile-acid-metabolizing enzymes
        and their small-molecule substrates, curated from public databases
        (UniProt, Rhea, ChEBI, RefSeq) and hand-audited literature recruitments.
      </p>

      <div className="mb-10 grid grid-cols-2 gap-3 sm:grid-cols-4">
        {[
          ["Proteins", S.n_proteins.toLocaleString()],
          ["Steroids", S.n_molecules.toLocaleString()],
          ["Nat + Syn", S.n_natsyn.toLocaleString()],
          ["Clusters", S.n_clusters.toString()],
        ].map(([label, val]) => (
          <div
            key={label}
            className="rounded-xl border border-slate-200 bg-white p-4 shadow-sm"
          >
            <div className="text-2xl font-bold text-slate-900">{val}</div>
            <div className="text-xs uppercase tracking-wide text-slate-500">
              {label}
            </div>
          </div>
        ))}
      </div>

      <h2 className="mb-4 text-2xl font-semibold text-slate-900">Pick a view</h2>
      <div className="grid gap-4 sm:grid-cols-3">
        <ViewCard
          href="/explore/protein"
          title="Protein centric"
          desc="Each point is a steroid-binding protein. Pick a cluster to see the family and its substrates."
          accent="from-sky-500 to-cyan-400"
        />
        <ViewCard
          href="/explore/molecule"
          title="Steroid centric"
          desc="Each point is a steroid molecule. Focus one to see the proteins that act on it."
          accent="from-amber-500 to-orange-400"
        />
        <ViewCard
          href="/explore/natsyn"
          title="Natural + synthetic"
          desc="Compare 2,889 natural and synthetic steroid entries in one shared chemical space."
          accent="from-emerald-500 to-lime-400"
        />
      </div>

      <footer className="mt-16 text-center text-xs text-slate-500">
        <span className="inline-flex items-center gap-1">
          <span className="text-amber-500">★</span>
          {S.n_new} newly-recruited from 2024–2026 literature.
        </span>{" "}
        · Built with Next.js + deck.gl.
      </footer>
    </main>
  );
}

function ViewCard({
  href,
  title,
  desc,
  accent,
}: {
  href: string;
  title: string;
  desc: string;
  accent: string;
}) {
  return (
    <Link
      href={href}
      className="group relative overflow-hidden rounded-2xl border border-slate-200 bg-white p-6 shadow-sm transition hover:-translate-y-0.5 hover:border-slate-300 hover:shadow-lg"
    >
      <div
        className={`absolute inset-x-0 top-0 h-1 bg-gradient-to-r ${accent}`}
      />
      <div className="text-lg font-semibold text-slate-900">{title}</div>
      <p className="mt-2 text-sm text-slate-600">{desc}</p>
      <div className="mt-4 text-sm font-medium text-sky-600 group-hover:underline">
        Explore →
      </div>
    </Link>
  );
}
