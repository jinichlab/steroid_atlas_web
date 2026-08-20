import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "Steroid Atlas",
  description:
    "Interactive UMAP atlas of steroid-metabolizing enzymes and their small-molecule substrates.",
};

export default function RootLayout({
  children,
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body className="min-h-screen bg-slate-50 text-slate-900 antialiased">
        {children}
      </body>
    </html>
  );
}
