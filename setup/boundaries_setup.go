package main

import (
	"context"
	"sync"
	"io"
	"strings"
	"fmt"
	"os/exec"
	"os"
	"encoding/json"
	"flag"
)

type Region struct {
	Chr string
	Start int
}

type Args struct {
	PairsPath string
	ChrSizesPath string
	Outpre string
	Binsize int
	Regions []Region
}

type Flags struct {
	Threads int
	SkipCload bool
	SkipBalance bool
	SkipZoom bool
	SkipBoundaries bool
	SkipPlot bool
}

func GetFlags() (Flags, error) {
	var f Flags

	flag.IntVar(&f.Threads, "t", 1, "Threads to use")
	flag.BoolVar(&f.SkipCload, "nocload", false, "Skip cload step")
	flag.BoolVar(&f.SkipBalance, "nobalance", false, "Skip balance step")
	flag.BoolVar(&f.SkipZoom, "nozoom", false, "Skip zoom step")
	flag.BoolVar(&f.SkipBoundaries, "noboundaries", false, "Skip boundary calling step")
	flag.BoolVar(&f.SkipPlot, "noplot", false, "Skip plotting step")

	flag.Parse()
	return f, nil
}

func GetArgs(r io.Reader) ([]Args, error) {
	var a Args
	var out []Args
	dec := json.NewDecoder(r)
	for e := dec.Decode(&a); e != io.EOF; e = dec.Decode(&a) {
		if e != nil {
			return nil, fmt.Errorf("GetArgs: %w", e)
		}
		out = append(out, a)
	}
	return out, nil
}

func SetCmdIo(cmd *exec.Cmd) {
	cmd.Stdin = os.Stdin
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
}

func Cload(ctx context.Context, a Args, threads int) error {
	cmd := exec.CommandContext(
		ctx,
		"cooler", "cload", "pairs",
		"-c1", "2",
		"-p1", "3",
		"-c2", "4",
		"-p2", "5",
		fmt.Sprintf("%s:%d", a.ChrSizesPath, a.Binsize),
		a.PairsPath,
		fmt.Sprintf("%s.%d.cool", a.Outpre, a.Binsize),
	)
	SetCmdIo(cmd)

	e := cmd.Run()
	if e != nil { return fmt.Errorf("Cload: %w", e) }
	return nil
}

func Balance(ctx context.Context, a Args, threads int) error {
	cmd := exec.CommandContext(
		ctx,
		"cooler", "balance",
		"-p", fmt.Sprintf("%d", threads),
		fmt.Sprintf("%s.%d.cool", a.Outpre, a.Binsize),
	)
	SetCmdIo(cmd)

	e := cmd.Run()
	if e != nil { return fmt.Errorf("Balance: %w", e) }
	return nil
}

func Zoom(ctx context.Context, a Args, threads int) error {
	cmd := exec.CommandContext(
		ctx,
		"cooler", "zoomify",
		"-p", fmt.Sprintf("%d", threads),
		"-r", "1000,10000,100000,1000000,10000000,100000000",
		fmt.Sprintf("%s.%d.cool", a.Outpre, a.Binsize),
	)
	SetCmdIo(cmd)

	e := cmd.Run()
	if e != nil { return fmt.Errorf("Zoom: %w", e) }
	return nil
}

func Boundaries(ctx context.Context, a Args, threads int) error {
	mcool := strings.ReplaceAll(
		fmt.Sprintf("%s.%d.cool", a.Outpre, a.Binsize),
		".cool", ".mcool",
	)

	cmdargs := []string {
		mcool,
		fmt.Sprintf("%s.boundaries", mcool),
	}
	for _, region := range a.Regions {
		cmdargs = append(cmdargs, fmt.Sprintf("%s:%d", region.Chr, region.Start))
	}

	cmd := exec.CommandContext(ctx, "boundaries_flexible.py", cmdargs...)
	SetCmdIo(cmd)

	e := cmd.Run()
	if e != nil { return fmt.Errorf("Boundaries: %w", e) }
	return nil
}

func Plot(ctx context.Context, a Args, threads int) error {
	return nil
}

func RunArg(ctx context.Context, a Args, f Flags) error {
	h := func(e error) error { return fmt.Errorf("RunArg: %w", e) }

	if !f.SkipCload {
		e := Cload(ctx, a, f.Threads)
		if e != nil { return h(e) }
	}

	if !f.SkipBalance {
		e := Balance(ctx, a, f.Threads)
		if e != nil { return h(e) }
	}

	if !f.SkipZoom {
		e := Zoom(ctx, a, f.Threads)
		if e != nil { return h(e) }
	}

	if !f.SkipBoundaries {
		e := Boundaries(ctx, a, f.Threads)
		if e != nil { return h(e) }
	}

	if !f.SkipPlot {
		e := Plot(ctx, a, f.Threads)
		if e != nil { return h(e) }
	}

	return nil
}

type Prunner struct {
	Jobs chan func(context.Context)
	Wg sync.WaitGroup
	Ctx context.Context
	Cancel context.CancelFunc
	CloseOnce sync.Once
}

func NewPrunner(threads int, ctx context.Context) (*Prunner) {
	p := new(Prunner)
	p.Ctx, p.Cancel = context.WithCancel(ctx)
	p.Jobs = make(chan func(context.Context))

	for i := 0; i < threads; i++ {
		go func() {
			defer func() {
				if r := recover(); r != nil {
					fmt.Println("recovering")
					p.Close()
					p.Cancel()
					fmt.Println("canceled")
					p.Wg.Done()
				}
			}()
			for job := range p.Jobs {
				job(p.Ctx)
				p.Wg.Done()
			}
		}()
	}
	return p
}

func (p *Prunner) Run(f func(context.Context)) {
	p.Wg.Add(1)
	p.Jobs <- f
}

func (p *Prunner) Wait() {
	p.Wg.Wait()
}

func (p *Prunner) Close() {
	p.CloseOnce.Do(func() {
		close(p.Jobs)
	})
}

func main() {
	flags, e := GetFlags()
	if e != nil { panic(e) }

	args, e := GetArgs(os.Stdin)
	if e != nil { panic(e) }

	p := NewPrunner(flags.Threads, context.Background())
	for _, arg := range args {
		arg := arg
		p.Run(func(ctx context.Context) {
			e := RunArg(ctx, arg, flags)
			if e != nil { panic(e) }
		})
		fmt.Println(arg)
	}
	p.Close()
	p.Wait()
}
